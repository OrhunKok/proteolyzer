import pandas as pd
import numpy as np
import os
from typing import Literal, Optional
import re
from proteoboost.utils.logging import MetaLogging

CELLEONE_MAPPING = {
                    'mTRAQ' : {
                                2 : {'G1' : '0', 'G2' : '8', 'G3' : '0', 'G4' : '8',
                                    'G5' : '0', 'G6' : '8', 'G7' : '0', 'G8' : '8'},

                                3 : {'G1' : '0', 'G2' : '4', 'G3' : '8', 'G4' : '0',
                                     'G5' : '4', 'G6' : '8', 'G7' : '0', 'G8' : '4',
                                     'G9' : '8', 'G10' : '0', 'G11' : '4', 'G12' : '8'}
                            },

                    'TMT' : {
                                14 : {'G1' : 'TMTpro-128C', 'G2' : 'TMTpro-129N', 'G3' : 'TMTpro-129C', 'G4' : 'TMTpro-130N',
                                    'G5' : 'TMTpro-130C', 'G6' : 'TMTpro-131N', 'G7' : 'TMTpro-131C', 'G8' : 'TMTpro-132N',
                                    'G9' : 'TMTpro-132C', 'G10' : 'TMTpro-133N', 'G11' : 'TMTpro-133C', 'G12' : 'TMTpro-134N',
                                    'G13' : 'TMTpro-134C', 'G14' : 'TMTpro-135N'}
                            }
                    }

class CoordinatesMapping(metaclass = MetaLogging):

    def __init__(
                self, 
                root_dir : str,
                label_type: Optional[Literal['mTRAQ', 'TMT']] = None,
                plex : Optional[int] = None
                ):
        
        pd.set_option('future.no_silent_downcasting', True)

        self.root_dir = root_dir
        self.label_type = label_type
        self.plex = plex
        self.file_paths = self._output_file_paths()
        self.parsed_data, self.parsed_stats = self._files_parse(self.file_paths)
        self.parsed_data = self._data_process(self.parsed_data)
        self.parsed_stats = self._stats_process(self.parsed_stats)

    
    def _output_file_paths(
                          self
                          ):
    
        geoprops_files, dispense_files, label_files, pickup_files = [], [], [], []
        for dirpath, dirnames, filenames in os.walk(self.root_dir):
            for filename in filenames:
        
                ### Geoprops Files ###
                if filename.endswith('.xls') and 'isolated' in filename and 'Reordered' in filename:
                    geoprops_files.append(os.path.join(dirpath, filename))     

                ### Dispense Files ###
                if filename.endswith('.log') and 'label' not in filename and 'pickup' not in filename:
                    dispense_files.append(os.path.join(dirpath, filename))  
                
                ### Label Files ###
                if filename.endswith('.log') and 'label' in filename:
                    label_files.append(os.path.join(dirpath, filename))  
                
                ### Pickup Files ###
                if filename.endswith('.log') and 'pickup' in filename:
                    pickup_files.append(os.path.join(dirpath, filename))
        
        cellenone_files = {'Geoprops' : geoprops_files, 'Dispense' : dispense_files, 'Label' : label_files, 'Pickup' : pickup_files}

        for key, value in cellenone_files.items():
            if len(value) == 0:
                print(f'No {key} files found in the directory: {self.root_dir}')
            else:
                print(f'Found {len(value)} {key} files.')

        cellenone_files = {key: value for key, value in cellenone_files.items() if len(value) > 0}

        return cellenone_files

    def _files_parse(
                    self,
                    file_paths : dict
                    ):

        parsed_data = {}
        parsed_stats = {}

        for key, value in file_paths.items():
            
            if key == 'Geoprops':
                df = self.xls_parse(value)
            elif key == 'Dispense':
                df, stats = self.log_parse(key, value)
            elif key == 'Label':
                df, stats = self.log_parse(key, value)
            elif key == 'Pickup':
                df, stats = self.log_parse(key, value)

            parsed_data[key] = df
            
            if key != 'Geoprops':
                parsed_stats[key] = stats

        return parsed_data, parsed_stats
    
    def _data_process(
                      self,
                      parsed_data : dict
                      ):

        for key, df in parsed_data.items():

            if key == 'Geoprops':
                df.columns = [re.sub(r'\s+', '', col) for col in df.columns]
                df['Timestamp'] = pd.to_datetime(df['Date'] + ' ' + df['Time'], format = 'mixed')

            if (key == 'Label') and (self.label_type is not None):
                df = self.label_well_plex(df)

            if key not in ['Geoprops', 'Dispense']:
                df['Timestamp'] = pd.to_datetime(df['Timestamp'], format = '%d.%m.%y-%H:%M:%S.%f', errors = 'coerce').dt.floor('s')
                df['Target'] = df['Target'].astype(int)
                df['Field'] = df['Field'].astype(int)
                df['XPos'] = df['XPos'].astype(int) + 1
                df['YPos'] = df['YPos'].astype(int) + 1


            if key == 'Pickup':

                # cellenOne only records the information for one nozzle. Other nozzles have to manually added as such.
                nozzle_well_mapping = {'A' : 'C', 'B' : 'D', 'E' : 'G', 'F' : 'H', 'I' : 'K', 'J' : 'L', 'M' : 'O', 'N' : 'P'}
                other_nozzle = df.copy()
                other_nozzle['Nozzle'] = 3
                other_nozzle['XPos'] = other_nozzle['XPos'] + 36
                other_nozzle['Well'] = other_nozzle['Well'].apply(lambda x: nozzle_well_mapping.get(x[0]) + x[1:])

                df = pd.concat([df, other_nozzle], axis = 0)
                df = df.reset_index(drop = True)

            parsed_data[key] = df

        return parsed_data


    def _metadata_validate(
                            self,
                            metadata : pd.DataFrame
                          ):
            
        if all(col in metadata.columns for col in ['Plate.Pickup', 'Well.Pickup']) and (self.plex is None):

            well_clash = metadata.groupby(['Plate.Pickup', 'Well.Pickup']).size() != 1
            metadata['Well.Clash'] = metadata.set_index(['Plate.Pickup', 'Well.Pickup']).index.map(well_clash).fillna(False)

            if well_clash.sum() > 0:
                self.logger.warning('Well clashes detected! Check Well.Clash column to view clashes.')


        if all(col in metadata.columns for col in ['Plate.Pickup', 'Well.Pickup', 'Plex.Label']):

            # metadata.groupby(['Plate.Pickup', 'Well.Pickup'])['Plex.Label'].apply(lambda x: x.duplicated(keep = False))
            labels_clash = metadata.groupby(['Plate.Pickup', 'Well.Pickup'])['Plex.Label'].nunique() != self.plex
            metadata['Label.Clash'] = metadata.set_index(['Plate.Pickup', 'Well.Pickup']).index.map(labels_clash).fillna(False)

            if labels_clash.sum() > 0:
                self.logger.warning('Label clashes detected! Check Label.Clash column to view clashes.')

        return metadata
    

    def map_data(self) -> pd.DataFrame:

        metadata = None

        if 'Pickup' in self.parsed_data and 'Geoprops' in self.parsed_data:
            pickup_df = self.parsed_data['Pickup']
            geo_df = self.parsed_data['Geoprops']

            pickup_df['MappedGeo'] = self._map_coords(geo_df, pickup_df)
            exploded = pickup_df.explode('MappedGeo')
            well_mapping = dict(zip(exploded['MappedGeo'], exploded['Well']))
            plate_mapping = dict(zip(exploded['MappedGeo'], exploded['Plate']))

            geo_df['Well.Pickup'] = geo_df.index.map(well_mapping)
            geo_df['Plate.Pickup'] = geo_df.index.map(plate_mapping)

            self.parsed_data['Geoprops'] = geo_df
            self.parsed_data['Pickup'] = pickup_df
            metadata = geo_df.copy()  

        if 'Geoprops' in self.parsed_data and 'Label' in self.parsed_data:
            geoprop_df = self.parsed_data['Geoprops'].copy()
            label_df = self.parsed_data['Label'].copy()

            geoprop_df = geoprop_df.rename(columns={
                col: f"{col}.Geoprops" for col in geoprop_df.columns
                if col not in ['XPos', 'YPos', 'Target', 'Field', 'Well.Pickup', 'Plate.Pickup']
            })
            label_df = label_df.rename(columns={
                col: f"{col}.Label" for col in label_df.columns
                if col not in ['XPos', 'YPos', 'Target', 'Field']
            })

            metadata = pd.merge(geoprop_df, label_df, on=['XPos', 'YPos', 'Target', 'Field'], how='inner')
        
        if metadata is not None:
            self._metadata_validate(metadata)
            return metadata
        else:
            self.logger.warning('Missing required files to do metadata mapping. Ensure all required inputs are present.')

        

    def _stats_process(
                      self,
                      parsed_stats : dict
                      ):

        for key, df in parsed_stats.items():
            
            df['Timestamp'] = df['Timestamp'].dt.floor('s')

            parsed_stats[key] = df

        return parsed_stats


    def map_stats(
                 self,
                 ):

        metastats = []
        for key, df in self.parsed_stats.items():
            
            df['Type'] = key

            metastats.append(df)

        metastats = pd.concat(metastats)

        return metastats

    
    def xls_parse(
                 self,
                 file_paths : list
                 ):
        
        dispense_dfs = []
        for file in file_paths:
            df = pd.read_csv(file, sep = '\t')
            df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

            # Replace with something smarter later
            sample_name = os.path.basename(os.path.dirname(file))
            sample_name = re.sub('.Run', '', sample_name)
            df['SampleName'] = np.repeat(sample_name, len(df))
            df = df.drop(['Plate', 'Well', 'ImageFile', 'Background'], axis = 1)

            dispense_dfs.append(df)

        return pd.concat(dispense_dfs, axis = 0).reset_index(drop = True)
    
    
    def log_parse(
                 self,
                 key : str,
                 file_paths : list
                 ):
        
        dfs = []
        stats_dfs = []
        for file in file_paths:
            with open(file, encoding = 'latin1') as file:
                df = file.readlines() 
                df = [re.sub('\n', '', ele) for ele in df]

                # Organize log file into dataframe
                df = pd.DataFrame(df)
                df = df[0].str.split('\t', expand = True)
                df = df.replace('', np.nan)
                df = df.dropna(how = 'all', axis = 1)

                temp_stats = df[(df.apply(lambda row: row.str.contains('Humidity')).any(axis=1)) & (df.apply(lambda row: row.str.contains('Temperature')).any(axis=1))]
                temp_stats = temp_stats.dropna(how = 'all', axis = 1)
                time_col = temp_stats.apply(lambda x: pd.to_datetime(x, format='%d.%m.%y-%H:%M:%S.%f', errors = 'coerce')).dropna(how = 'all', axis = 1)
                val_cols = temp_stats.apply(lambda x: pd.to_numeric(x, errors = 'coerce')).dropna(how = 'all', axis = 1)
                temp_stats = pd.concat([time_col, val_cols], axis = 1).reset_index(drop = True)
                temp_stats.columns = ['Timestamp', 'Humidity', 'Temperature', 'Dew Point', 'Adj. Temp', 'Bath Temp']
                temp_stats['Timestamp'] = temp_stats['Timestamp'].dt.floor('s')
                stats_dfs.append(temp_stats)

                # Subset the log messages (relevant for label and pickup)
                df = df[df.apply(lambda row: row.str.contains('drops')).any(axis=1)]
                df = df.dropna(how = 'all', axis = 1)

                if len(df) == 0:
                    # Can exit from dispense log here
                    if key == 'Dispense':
                        dfs.append(pd.DataFrame())
                    else:
                        # In case log from aborted process is included in the output
                        print(f'Skipping: {file.name}. Likely an aborted process during sample-prep.')
                
                else:
                    df.columns = ['Timestamp', 'Plate', 'PlatePos', 'Nozzle', 'Well', 'Target', 'Level', 'Field', 'Drops', 'XPos', 'YPos']
                    df = df.reset_index(drop = True)

                    if key == 'Label':
                        df = df.drop(['PlatePos'], axis = 1)

                    # Quick fix for plate number for now
                    if key == 'Pickup':
                        pickup_name = os.path.basename(file.name)
                        plate_num = re.findall(r'pickup_(.*?)_', pickup_name)
                        df['Plate'] = int(plate_num[0])

                    dfs.append(df)

        dfs = pd.concat(dfs, axis = 0)
        stats_dfs = pd.concat(stats_dfs, axis = 0)
        return dfs, stats_dfs

    def label_well_plex(
                        self,
                        label_df : pd.DataFrame
                        ):
        
        label_mapping = CELLEONE_MAPPING[self.label_type][self.plex]
        label_df['Plex'] = label_df['Well'].map(label_mapping)

        return label_df
    
    def _map_coords(
                    self, 
                    geo_df,
                    map_df, 
                    coord_cols=['XPos', 'YPos'], 
                    group_cols=['Target', 'Field']
                    ):

        grouped_geo = geo_df.groupby(group_cols)
        grouped_map = map_df.groupby(group_cols)
        
        results = []
        for group_key, geo_group in grouped_geo:
            if group_key not in grouped_map.groups:
                print(f'{group_key} not in pickup groups, skipping.')
                continue
            
            map_group = grouped_map.get_group(group_key)

            geo_coords = np.stack(geo_group[coord_cols].values)
            map_coords = np.stack(map_group[coord_cols].values)

            distances = np.linalg.norm(map_coords[:, np.newaxis, :] - geo_coords[np.newaxis, :, :], axis=2)
            sorted_indices = np.argsort(distances, axis=1)[:, :self.plex]
            closest_points = geo_group.index.values[sorted_indices]

            for i, map_idx in enumerate(map_group.index):
                results.append({
                    'closest_points': closest_points[i],
                })

        return pd.DataFrame(results)