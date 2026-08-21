CELLEONE_MAPPING = {
    "mTRAQ": {
        2: {
            "G1": "0",
            "G2": "8",
            "G3": "0",
            "G4": "8",
            "G5": "0",
            "G6": "8",
            "G7": "0",
            "G8": "8",
        },
        3: {
            "G1": "0",
            "G2": "4",
            "G3": "8",
            "G4": "0",
            "G5": "4",
            "G6": "8",
            "G7": "0",
            "G8": "4",
            "G9": "8",
            "G10": "0",
            "G11": "4",
            "G12": "8",
        },
    },
    "TMT": {
        14: {
            "G1": "TMTpro-128C",
            "G2": "TMTpro-129N",
            "G3": "TMTpro-129C",
            "G4": "TMTpro-130N",
            "G5": "TMTpro-130C",
            "G6": "TMTpro-131N",
            "G7": "TMTpro-131C",
            "G8": "TMTpro-132N",
            "G9": "TMTpro-132C",
            "G10": "TMTpro-133N",
            "G11": "TMTpro-133C",
            "G12": "TMTpro-134N",
            "G13": "TMTpro-134C",
            "G14": "TMTpro-135N",
        }
    },
}


# cellenOne only records the information for one nozzle. Other nozzles have to be manually added as such.
NOZZLE_WELL_MAPPING = {
    "A": "C",
    "B": "D",
    "E": "G",
    "F": "H",
    "I": "K",
    "J": "L",
    "M": "O",
    "N": "P",
}

PICKUP_NOZZLE_ID = 3
PICKUP_NOZZLE_XPOS_OFFSET = 36


#: Geoprops column naming the imaging channel a measurement came from. The
#: export writes one row per (cell, channel), e.g. Transmission and Green.
CHANNEL_COL = "Teg"

#: Channel whose measurements keep the unsuffixed column names when several are
#: present; the others become e.g. "Diameter.Green". Transmission is the
#: detection channel, so it is the one that always carries the geometry.
PRIMARY_CHANNEL = "Transmission"


TEMP_STATS_COLS = [
    "Timestamp",
    "Humidity",
    "Temperature",
    "Dew Point",
    "Adj. Temp",
    "Bath Temp",
]

DROPLET_COLS = [
    "Timestamp",
    "Plate",
    "PlatePos",
    "Nozzle",
    "Well",
    "Target",
    "Level",
    "Field",
    "Drops",
    "XPos",
    "YPos",
]


MERGE_COLS = ["XPos", "YPos", "Target", "Field"]
