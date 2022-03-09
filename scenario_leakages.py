# -*- coding: utf-8 -*-
"""
WDN Dataset Generator Leakage scenarios
Copyright: (C) 2022, KIOS Research Center of Excellence
"""
import pandas as pd
from numpy import isnan, arange
import wntr
import pickle
import os
import sys
import yaml
import shutil
import time
from math import sqrt
import os
import logging

# Read input arguments from yalm file
try:
    with open(os.path.join(os.getcwd(), 'dataset_configuration.yalm'), 'r') as f:
        leak_pipes = yaml.full_load(f.read())
except:
    print('"dataset_configuration" file not found.')
    logging.info('"dataset_configuration" file not found.')
    os.startfile(logfilename)
    sys.exit(1)

def get_values(leak_pipes, field):
    values = []
    [values.append(str(sens)) for sens in leak_pipes[field] if sens is not None]
    return values


leakages = pd.read_excel('create_leakage_scenarios.xlsx', engine='openpyxl', converters={'scenario': int, 'linkid': str})

start_time = leak_pipes['times']['StartTime']
end_time = leak_pipes['times']['EndTime']
inp_file = leak_pipes['Network']['filename']
results_folder = f'{os.getcwd()}\\LeakageScenarios\\'

pressure_sensors = get_values(leak_pipes, 'pressure_sensors')
amrs = get_values(leak_pipes, 'amrs')
flow_sensors = get_values(leak_pipes, 'flow_sensors')
level_sensors = get_values(leak_pipes, 'level_sensors')

# Check if sensor exists
logfilename = "dataset_generator.log"
errcode = False
logging.basicConfig(filename=logfilename, level=logging.INFO, filemode="w")
print(f'Run input file: "{inp_file}"')
logging.info(f'Run input file: "{inp_file}"')
logging.info('Start leakage scenarios dataset generator.')
logging.info('Check configuration yalm file.')

# demand-driven (DD) or pressure dependent demand (PDD)
Mode_Simulation = 'PDD'  # 'PDD'#'PDD'


class LeakDatasetCreator:
    def __init__(self, scenario):

        # Create Results folder
        if scenario == 1:
            self.create_folder(results_folder)
        self.unc_range = arange(0, 0.25, 0.05)

        self.scenario = scenario
        # Load EPANET network file
        self.wn = wntr.network.WaterNetworkModel(inp_file)

        for name, node in self.wn.junctions():
            node.required_pressure = 25

        self.inp = os.path.basename(self.wn.name)[0:-4]

        # Get the name of input file
        self.net_name = f'{results_folder}{self.inp}'

        # Get time step
        self.time_step = round(self.wn.options.time.hydraulic_timestep)
        # Create time_stamp
        try:
            self.time_stamp = pd.date_range(start_time, end_time, freq=str(self.time_step / 60) + "min")
        except:
            print('Please check you time step in network file.')
            logging.info('Please check you time step in network file.')
            os.startfile(logfilename)
            sys.exit(1)

        # Simulation duration in steps
        self.wn.options.time.duration = (len(self.time_stamp) - 1) * 300  # 5min step
        self.TIMESTEPS = int(self.wn.options.time.duration / self.wn.options.time.hydraulic_timestep)

    def create_csv_file(self, values, time_stamp, columnname, pathname):

        file = pd.DataFrame(values)
        file['time_stamp'] = time_stamp
        file = file.set_index(['time_stamp'])
        file.columns.values[0] = columnname
        file.to_csv(pathname)
        del file, time_stamp, values

    def create_folder(self, _path_):

        try:
            if os.path.exists(_path_):
                shutil.rmtree(_path_)
            os.makedirs(_path_)
        except Exception as error:
            pass

    def dataset_generator(self):
        # Path of EPANET Input File
        print(f"Generating dataset...")
        logging.info(f"Generating dataset...")

        # Initialize parameters for the leak
        leak_node = {}
        leak_diameter = {}
        leak_area = {}
        leak_type = {}
        leak_starts = {}
        leak_ends = {}
        leak_peak_time = {}
        leak_param = {}
        leak_i = 0

        number_of_leaks = list(set(leakages['scenario'])).count(self.scenario)
        scenario_rows = leakages.iloc[:, 0] == scenario
        scenario_rows = [i for i, x in enumerate(scenario_rows) if x]
        for leakn in [leakages.iloc[scenario_rows]]:
            # Split pipe and add a leak node
            # leakages: pipeID, startTime, endTime, leakDiameter, leakType (abrupt, incipient)
            # Start time of leak
            #leakages.iloc[:, 0] == scenario
            lind = leakn.index.values[0]
            if leakn['scenario'][lind] != self.scenario:
                continue
            ST = self.time_stamp.get_loc(leakn['starttime'][lind])

            # End Time of leak
            ET = self.time_stamp.get_loc(leakn['endtime'][lind])

            # Get leak type
            leak_type[leak_i] = leakn['profile'][lind]

            # Split pipe to add a leak
            pipe_id = self.wn.get_link(str(leakn['linkid'][lind]))
            node_leak = f'{pipe_id}_leaknode'
            self.wn = wntr.morph.split_pipe(self.wn, pipe_id, f'{pipe_id}_Bleak', node_leak)
            leak_node[leak_i] = self.wn.get_node(self.wn.node_name_list[self.wn.node_name_list.index(node_leak)])

            if 'incipient' in leakn['profile'][lind]:
                # END TIME
                ET = ET + 1
                PT = self.time_stamp.get_loc(leakn['peaktime'][lind])+1

                # Leak diameter as max magnitude for incipient
                nominal_pres = 100
                leak_diameter[leak_i] = float(leakn['leakdiameter'][lind])
                leak_area[leak_i] = 3.14159 * (leak_diameter[leak_i] / 2) ** 2

                # incipient
                leak_param[leak_i] = 'demand'
                increment_leak_diameter = leak_diameter[leak_i] / (PT - ST)
                increment_leak_diameter = arange(increment_leak_diameter, leak_diameter[leak_i], increment_leak_diameter)
                increment_leak_area = 0.75 * sqrt(2 / 1000) * 990.27 * 3.14159 * (increment_leak_diameter/2)**2
                leak_magnitude = 0.75 * sqrt(2 / 1000) * 990.27 * leak_area[leak_i]
                pattern_array = [0] * (ST) + increment_leak_area.tolist() + [leak_magnitude] * (ET - PT + 1) + [0] * (self.TIMESTEPS - ET)

                # basedemand
                leak_node[leak_i].demand_timeseries_list[0]._base = 1
                pattern_name = f'{str(leak_node[leak_i])}'
                self.wn.add_pattern(pattern_name, pattern_array)
                leak_node[leak_i].demand_timeseries_list[0].pattern_name = pattern_name
                leak_node[leak_i].required_pressure = nominal_pres
                leak_node[leak_i].minimum_pressure = 0

                # save times of leak
                leak_starts[leak_i] = self.time_stamp[ST]
                leak_starts[leak_i] = leak_starts[leak_i]._date_repr + ' ' + leak_starts[leak_i]._time_repr
                leak_ends[leak_i] = self.time_stamp[ET - 1]
                leak_ends[leak_i] = leak_ends[leak_i]._date_repr + ' ' + leak_ends[leak_i]._time_repr
                leak_peak_time[leak_i] = self.time_stamp[PT - 1]._date_repr + ' ' + self.time_stamp[PT - 1]._time_repr

            else:
                leak_param[leak_i] = 'leak_demand'
                PT = ST
                leak_diameter[leak_i] = float(leakn['leakdiameter'][lind])
                leak_area[leak_i] = 3.14159 * (leak_diameter[leak_i] / 2) ** 2

                leak_node[leak_i]._leak_end_control_name = str(leak_i) + 'end'
                leak_node[leak_i]._leak_start_control_name = str(leak_i) + 'start'

                leak_node[leak_i].add_leak(self.wn, discharge_coeff=0.75,
                                           area=leak_area[leak_i],
                                           start_time=ST * self.time_step,
                                           end_time=(ET+1) * self.time_step)

                leak_starts[leak_i] = self.time_stamp[ST]
                leak_starts[leak_i] = leak_starts[leak_i]._date_repr + ' ' + leak_starts[leak_i]._time_repr
                leak_ends[leak_i] = self.time_stamp[ET]
                leak_ends[leak_i] = leak_ends[leak_i]._date_repr + ' ' + leak_ends[leak_i]._time_repr
                leak_peak_time[leak_i] = self.time_stamp[PT]._date_repr + ' ' + self.time_stamp[PT]._time_repr
            leak_i += 1

        # Save/Write input file with new settings
        if number_of_leaks:
            leakages_folder = f'{results_folder}scenario{str(self.scenario)}\\Leakages'
            self.create_folder(leakages_folder)

        # Save the water network model to a file before using it in a simulation
        with open('self.wn.pickle_leak', 'wb') as f:
            pickle.dump(self.wn, f)

        # Run wntr simulator
        self.wn.options.hydraulic.demand_model = Mode_Simulation
        sim = wntr.sim.WNTRSimulator(self.wn)
        self.results = sim.run_sim()

        if self.results.node["pressure"].empty:
            print("Negative pressures.")
            logging.info("Negative pressures.")
            return -1

        if self.results:
            # Create CSV files
            for leak_i in range(0, len(leak_node)):

                if 'leaknode' in str(leak_node[leak_i]):
                    NODEID = str(leak_node[leak_i]).split('_')[0]
                totals_info = {'Description': ['Leak Pipe', 'Leak Area', 'Leak Diameter', 'Leak Type', 'Leak Start',
                                               'Leak End', 'Peak Time'],
                               'Value': [NODEID, str(leak_area[leak_i]), str(leak_diameter[leak_i]),
                                         leak_type[leak_i],
                                         str(leak_starts[leak_i]), str(leak_ends[leak_i]), str(leak_peak_time[leak_i])]}

                # Create leak XLS files
                decimal_size = 2

                leaks = self.results.node[leak_param[leak_i]][str(leak_node[leak_i])].values
                # Convert m^3/s (wntr default units) to m^3/h
                # https://wntr.readthedocs.io/en/latest/units.html#epanet-unit-conventions
                leaks = [elem * 3600 for elem in leaks]
                leaks = [round(elem, decimal_size) for elem in leaks]
                leaks = leaks[:len(self.time_stamp)]

                total_Leaks = {'Timestamp': self.time_stamp}
                total_Leaks[NODEID] = leaks
                #self.create_csv_file(leaks, self.time_stamp, 'Description', f'{results_folder}scenario{str(scenario)}\\Leak_{str(leak_node[leak_i])}_demand.csv')
                df1 = pd.DataFrame(totals_info)
                df2 = pd.DataFrame(total_Leaks)
                writer = pd.ExcelWriter(f'{leakages_folder}\\Leak_{NODEID}.xlsx', engine='xlsxwriter')
                df1.to_excel(writer, sheet_name='Info', index=False)
                df2.to_excel(writer, sheet_name='Demand (m3_h)', index=False)
                writer.save()

            # Create xlsx file with Measurements
            total_pressures = {'Timestamp': self.time_stamp}
            total_demands = {'Timestamp': self.time_stamp}
            total_flows = {'Timestamp': self.time_stamp}
            total_levels = {'Timestamp': self.time_stamp}
            for j in range(0, self.wn.num_nodes):
                node_id = self.wn.node_name_list[j]

                if node_id in pressure_sensors:
                    pres = self.results.node['pressure'][node_id]
                    pres = pres[:len(self.time_stamp)]
                    pres = [round(elem, decimal_size) for elem in pres]
                    total_pressures[node_id] = pres

                if node_id in amrs:
                    dem = self.results.node['demand'][node_id]
                    dem = dem[:len(self.time_stamp)]
                    dem = [elem * 3600 * 1000 for elem in dem] #CMH / L/s
                    dem = [round(elem, decimal_size) for elem in dem]
                    total_demands[node_id] = dem

                if node_id in level_sensors:
                    level_pres = self.results.node['pressure'][node_id]
                    level_pres = level_pres[:len(self.time_stamp)]
                    level_pres = [round(elem, decimal_size) for elem in level_pres]
                    total_levels[node_id] = level_pres

            for j in range(0, self.wn.num_links):
                link_id = self.wn.link_name_list[j]
                if link_id not in flow_sensors:
                    continue
                flows = self.results.link['flowrate'][link_id]
                flows = flows[:len(self.time_stamp)]
                flows = [elem * 3600 for elem in flows]
                flows = [round(elem, decimal_size) for elem in flows]
                total_flows[link_id] = flows

            # Create a Pandas dataframe from the data.
            df1 = pd.DataFrame(total_pressures)
            df2 = pd.DataFrame(total_demands)
            df3 = pd.DataFrame(total_flows)
            df4 = pd.DataFrame(total_levels)
            # Create a Pandas Excel writer using XlsxWriter as the engine.
            writer = pd.ExcelWriter(f'{results_folder}scenario{str(self.scenario)}\\Measurements.xlsx', engine='xlsxwriter')

            # Convert the dataframe to an XlsxWriter Excel object.
            # Pressures (m), Demands (m^3/h), Flows (m^3/h), Levels (m)
            df1.to_excel(writer, sheet_name='Pressures (m)', index=False)
            df2.to_excel(writer, sheet_name='Demands (L_h)', index=False)
            df3.to_excel(writer, sheet_name='Flows (m3_h)', index=False)
            df4.to_excel(writer, sheet_name='Levels (m)', index=False)

            # Close the Pandas Excel writer and output the Excel file.
            writer.save()

            try:
                os.remove('self.wn.pickle_leak')
            except:
                pass
        else:
            print('Results empty.')
            logging.info('Results empty.')
            return -1


if __name__ == '__main__':

    # Create tic / toc
    t = time.time()

    # Call leak dataset creator
    for scenario in leakages['scenario']:
        if isnan(scenario):
            continue
        L = LeakDatasetCreator(scenario)
        L.dataset_generator()

    print(f"Dataset completed.")
    logging.info(f"Dataset completed.")
    print(f'Total Elapsed time is {str(time.time() - t)} seconds.')
    logging.info(f'Total Elapsed time is {str(time.time() - t)} seconds.')
    os.startfile(logfilename)