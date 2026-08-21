"""Run the AAS stages in order, in the two phases the method allows.

The five stages are not one automated run: detection writes a FASTA that has
to be searched against the raw files by the search engine, and only then do
the ``*_val`` directories the validation phase reads exist. So there are two
phases with a manual step between them:

    pipeline = aas.Pipeline("params.yaml")

    pipeline.run_detection()    # preprocess -> translate -> detect
    #   ... search the raw files against <output>/<sample>_validation.fasta ...
    pipeline.run_validation()   # preprocess -> validate -> quantify

    pipeline.status()           # what has run, what can run now

What this buys over calling the stages directly is the ordering: the
preprocessor runs again in the second phase to convert the validation
searches, translation is skipped when its frames already exist, and phase two
refuses to start before the searches are there.
"""

from pathlib import Path

from proteolyzer.core.logging import Logged

from .detection import Detection
from .params import load_params
from .preprocessing import Preprocessor
from .quantification import Quantification
from .results import Results
from .translation import FrameTranslator
from .validation import Validation

#: Suffix cellenONE-independent validation search directories carry.
VALIDATION_SUFFIX = "_val"


class Pipeline(Logged):
    """The AAS stages for one parameter file, run in the right order."""

    def __init__(self, params, queue=None):
        self.raw_params = params
        self.params = load_params(params)
        self.queue = queue
        self.output_dir = Path(self.params["Utils"]["Output Folder"])
        self.data_dir = Path(self.params["Utils"]["Data Folder"])
        self.frames_dir = Path(self.params["Translation"]["Translated Frames Folder"])

    # ------------------------------------------------------------ inspection

    @property
    def results(self) -> Results:
        """What this run has produced so far."""
        return Results(self.output_dir)

    def frames_translated(self) -> bool:
        """Whether the six-frame translation has already been written."""
        return all(
            (self.frames_dir / f"frame_{frame}.p").exists() for frame in range(1, 7)
        )

    def validation_searches(self) -> list[Path]:
        """The validation search directories, once the manual search has run."""
        if not self.data_dir.is_dir():
            return []
        return sorted(
            path
            for path in self.data_dir.rglob(f"*{VALIDATION_SUFFIX}")
            if path.is_dir()
        )

    def status(self) -> dict:
        """What has run and what can run now."""
        searches = self.validation_searches()
        return {
            "frames_translated": self.frames_translated(),
            "validation_searches": [path.name for path in searches],
            "samples": self.results.samples,
            "can_run_detection": self.data_dir.is_dir(),
            "can_run_validation": bool(searches),
        }

    # ----------------------------------------------------------------- phases

    def run_detection(self, translate: bool | None = None) -> Results:
        """Phase one: preprocess, translate the genome, detect candidates.

        Parameters
        ----------
        translate
            Whether to run the six-frame translation. The default skips it when
            its frames are already on disk, since it is the slow step and its
            output only depends on the genome.
        """
        self.logger.info("Phase one: preprocessing, translation, detection")
        Preprocessor.MaxQuant(self.raw_params, self.queue).run()

        if translate is None:
            translate = not self.frames_translated()
        if translate:
            FrameTranslator(self.raw_params, self.queue).run()
        else:
            self.logger.info(
                f"Translated frames already in {self.frames_dir}; skipping "
                "translation. Pass translate=True to redo it."
            )

        Detection(self.raw_params, self.queue).run()

        self.logger.info(
            "Phase one done. Search the raw files against the validation FASTA "
            f"in {self.output_dir}, then call run_validation()."
        )
        return self.results

    def run_validation(self) -> Results:
        """Phase two: convert the validation searches, validate, quantify."""
        searches = self.validation_searches()
        if not searches:
            raise FileNotFoundError(
                f"No '*{VALIDATION_SUFFIX}' search directories under "
                f"{self.data_dir}. Phase two reads the search of the validation "
                "FASTA that run_detection() wrote; run that search first."
            )

        self.logger.info(
            f"Phase two: {len(searches)} validation search(es), then validation "
            "and quantification"
        )
        # Again, because the validation searches did not exist the first time.
        Preprocessor.MaxQuant(self.raw_params, self.queue).run()
        Validation(self.raw_params, self.queue).run()
        Quantification(self.raw_params, self.queue).run()

        self.logger.info("Phase two done.")
        return self.results
