# -*- coding: utf-8 -*-
"""
Multi WDN Dataset Generator Sensor Faults
Copyright: (C) 2022, KIOS Research Center of Excellence
"""
import pandas as pd
from numpy import exp, random, arange
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
        data_yalm = yaml.full_load(f.read())
except:
    print('"dataset_configuration" file not found.')
    logging.info('"dataset_configuration" file not found.')
    os.startfile(logfilename)
    sys.exit(1)

def get_values(data_yalm, field):
    values = []
    [values.append(str(sens)) for sens in data_yalm[field] if sens is not None]
    return values


sensor_faults = pd.read_excel('create_sensor_fault_scenarios.xlsx', engine='openpyxl')

start_time = data_yalm['times']['StartTime']
end_time = data_yalm['times']['EndTime']
inp_file = data_yalm['Network']['filename']
scenariofolder = 'SensorFaultScenarios'
results_folder = f'{os.getcwd()}\\{scenariofolder}\\'

pressure_sensors = get_values(data_yalm, 'pressure_sensors')
amrs = get_values(data_yalm, 'amrs')
flow_sensors = get_values(data_yalm, 'flow_sensors')
level_sensors = get_values(data_yalm, 'level_sensors')
sensortypes = {'level': 'pressure', 'flow': 'flowrate', 'pressure': 'pressure', 'amrs': 'demand'}
objecttype = {'level': 'node', 'flow': 'link', 'pressure': 'node', 'amrs': 'node'}
sheetnamefault = {'level': 'Levels (m)', 'flow': 'Flows (m3_h)', 'pressure': 'Pressures (m)', 'amrs': 'Demands (L_h)'}

# Check if sensor exists
logfilename = "dataset_generator.log"
errcode = False
logging.basicConfig(filename=logfilename, level=logging.INFO, filemode="w")
print(f'Run input file: "{inp_file}"')
logging.info(f'Run input file: "{inp_file}"')
logging.info('Start sensor fault scenarios dataset generator.')
logging.info('Check configuration yalm file.')

for sfault in sensor_faults.iterrows():
    id = str(sfault[1]['nodeid/linkid'])
    typefault = sfault[1]['sensortype']
    if typefault == 'pressure' and id not in pressure_sensors:
        errcode = True
        logging.error(f'Pressure sensor "{id}" does not exist at the location of the sensor fault!')

    if typefault == 'flow' and id not in flow_sensors:
        errcode = True
        logging.error(f'Flow sensor "{id}" does not exist at the location of the sensor fault!')

    if typefault == 'level' and id not in level_sensors:
        errcode = True
        logging.error(f'Level sensor "{id}" does not exist at the location of the sensor fault!')

    if typefault == 'amrs' and id not in amrs:
        errcode = True
        logging.error(f'Amrs sensor "{id}" does not exist at the location of the sensor fault!')

if errcode:
    print('Error: A sensor does not exist at the location of the sensor fault!')
    logging.info('Stop script generator.')
    os.startfile(logfilename)
    sys.exit(1)

# demand-driven (DD) or pressure dependent demand (PDD)
Mode_Simulation = 'PDD'  # 'PDD'#'PDD'


class DatasetCreator:
    def __init__(self):

        # Create Results folder
        self.create_folder(results_folder)

        self.scenario_num = 1
        self.unc_range = arange(0, 0.25, 0.05)

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
        self.runsimulator = True

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

    def sensorfaultmodels(self, y0, index_id, fstart, fend, ftype, fpar, a1, a2):
        # Sensor faults: https://github.com/eldemet/sensorfaultmodels/blob/main/sensorfaultmodels.m
        #                https://github.com/Mariosmsk/sensorfaultmodels/blob/main/sensorfaultmodels.py
        T1 = fstart[index_id]
        T2 = fend[index_id]
        a1 = a1[index_id]
        a2 = a2[index_id]
        ftype = ftype[index_id]
        fpar = float(fpar[index_id])

        y = []
        for k in range(0, len(y0)):
            y0k = y0[k]
            b1 = 0
            b2 = 0
            if k >= T1:
                b1 = 1 - exp(- a1 * (k - T1))

            if k >= T2:
                b2 = 1 - exp(- a2 * (k - T2))

            b = b1 - b2
            phi = 0

            if b > 0:
                if ftype == 'constant':
                    phi = fpar
                if ftype == 'drift':
                    phi = fpar * (k - T1)
                if ftype == 'normal':
                    phi = random.normal(0, fpar)
                if ftype == 'percentage':
                    phi = fpar * y0k
                if ftype == 'stuckzero':
                    phi = -y0k

            df = b * phi
            y0k = y0k + df
            y.append(y0k)
        return y

    def dataset_generator(self, scenario):
        # Path of EPANET Input File
        print(f"Generating dataset...")
        logging.info(f"Generating dataset...")

        # Save the water network model to a file before using it in a simulation
        with open('self.wn.pickle_fault', 'wb') as f:
            pickle.dump(self.wn, f)

        # Run wntr simulator
        if self.runsimulator:
            self.wn.options.hydraulic.demand_model = Mode_Simulation
            sim = wntr.sim.WNTRSimulator(self.wn)
            self.results = sim.run_sim()
            self.runsimulator = False

        # Faults
        fault_objid = {}
        fault_param = {}
        fault_type = {}
        fault_functionpar = {}
        fault_sens_start = {}
        fault_sens_end = {}
        fault_a1 = {}
        fault_a2 = {}
        fault_index_ids = {}
        fault_scenario_id = {}
        for sfault in sensor_faults.iterrows():

            # Get sensor fault
            scenario_id = sfault[1]['scenario']
            if scenario_id != scenario:
                continue
            fault_id = str(sfault[1]['nodeid/linkid'])
            sens_type = sfault[1]['sensortype']
            functiontype = sfault[1]['profile']
            functionpar = sfault[1]['magnitude']
            fault_start = sfault[1]['starttime']
            fault_end = sfault[1]['endtime']
            index_id = sens_type + fault_id
            fault_index_ids[index_id] = index_id

            # Start time of leak
            ST = self.time_stamp.get_loc(fault_start)

            # End Time of leak
            ET = self.time_stamp.get_loc(fault_end)

            # Split pipe to add a leak
            if objecttype[sens_type] == 'link':
                object_id = self.wn.get_link(fault_id)
            else:
                object_id = self.wn.get_node(fault_id)

            if len(sensor_faults) > 0:
                if 'constant' in functiontype:
                    fault_a1[index_id] = 0.5  # parameter in occurance evolution profile function
                    fault_a2[index_id] = 0.7  # parameter in dissapearance evolution profile function
                if 'drift' in functiontype:
                    fault_a1[index_id] = 999999
                    fault_a2[index_id] = 1
                if 'normal' in functiontype:
                    fault_a1[index_id] = 100
                    fault_a2[index_id] = 100
                if 'percentage' in functiontype:
                    fault_a1[index_id] = 999999
                    fault_a2[index_id] = .7
                if 'stuckzero' in functiontype:
                    fault_a1[index_id] = 999999
                    fault_a2[index_id] = .7

                fault_objid[index_id] = object_id.name
                fault_param[index_id] = sens_type
                fault_scenario_id[index_id] = scenario_id
                fault_type[index_id] = functiontype
                fault_functionpar[index_id] = functionpar
                fault_sens_start[index_id] = ST
                fault_sens_end[index_id] = ET

        # Save/Write input file with new settings
        if self.results.node["pressure"].empty:
            print("Negative pressures.")
            logging.info("Negative pressures.")
            return -1

        if self.results:
            # Create CSV files
            faults_folder = f'{results_folder}scenario{str(scenario)}\\WithoutSensorFaults'
            self.create_folder(faults_folder)
            for index_id in fault_index_ids.values():
                fault_id = fault_objid[index_id]
                T1 = fault_sens_start[index_id]
                T2 = fault_sens_end[index_id]
                ftype = fault_type[index_id]
                fpar = float(fault_functionpar[index_id])

                str_T1 = self.time_stamp[T1]
                str_T1 = str_T1._date_repr + ' ' + str_T1._time_repr

                str_T2 = self.time_stamp[T2]
                str_T2 = str_T2._date_repr + ' ' + str_T2._time_repr

                totals_info = {'Description': [f'{objecttype[fault_param[index_id]].upper()} ID', 'Sensor Type', 'Function type', 'Function par', 'Fault Start',
                                               'Fault End'],
                               'Value': [index_id, str(sensortypes[fault_param[index_id]]), str(ftype),
                                         str(fpar), str(str_T1), str(str_T2)]}
                # Create faults XLS files
                decimal_size = 2

                withoutfault = eval(f"self.results.{objecttype[fault_param[index_id]]}['{sensortypes[fault_param[index_id]]}']['{fault_id}'].values")
                withoutfault = withoutfault[:len(self.time_stamp)]
                if 'demand' in sensortypes[fault_param[index_id]]:
                    withoutfault = [elem * 3600 * 1000 for elem in withoutfault] #CMH / L/s
                elif 'flow' in sensortypes[fault_param[index_id]]:
                    withoutfault = [elem * 3600 for elem in withoutfault]

                # Convert m^3/s (wntr default units) to m^3/h
                # https://wntr.readthedocs.io/en/latest/units.html#epanet-unit-conventions
                withoutfault = [round(elem, decimal_size) for elem in withoutfault]
                withoutfault = withoutfault[:len(self.time_stamp)]

                total_Faults = {'Timestamp': self.time_stamp}
                total_Faults[fault_id] = withoutfault
                df1 = pd.DataFrame(totals_info)
                df2 = pd.DataFrame(total_Faults)
                writer = pd.ExcelWriter(f'{faults_folder}\\WithoutSensorFault_{index_id}_{fault_param[index_id]}.xlsx', engine='xlsxwriter')
                df1.to_excel(writer, sheet_name='Info', index=False)
                df2.to_excel(writer, sheet_name=f'{sheetnamefault[fault_param[index_id]]}', index=False)
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
                    if f'pressure{node_id}' in fault_objid.keys():
                        index_id = fault_index_ids[f'pressure{node_id}']
                        pres = self.sensorfaultmodels(pres.values, index_id, fault_sens_start, fault_sens_end,
                                                      fault_type, fault_functionpar, fault_a1, fault_a2)

                    pres = [round(elem, decimal_size) for elem in pres]
                    total_pressures[node_id] = pres

                if node_id in amrs:
                    dem = self.results.node['demand'][node_id]
                    dem = dem[:len(self.time_stamp)]
                    dem = [elem * 3600 * 1000 for elem in dem] #CMH / L/s
                    if f'amrs{node_id}' in fault_objid.keys():
                        index_id = fault_index_ids[f'amrs{node_id}']
                        dem = self.sensorfaultmodels(dem, index_id, fault_sens_start, fault_sens_end,
                                                      fault_type, fault_functionpar, fault_a1, fault_a2)
                    dem = [round(elem, decimal_size) for elem in dem]
                    total_demands[node_id] = dem

                if node_id in level_sensors:
                    level_pres = self.results.node['pressure'][node_id]
                    level_pres = level_pres[:len(self.time_stamp)]
                    if f'level{node_id}' in fault_objid.keys():
                        index_id = fault_index_ids[f'level{node_id}']
                        level_pres = self.sensorfaultmodels(level_pres.values, index_id, fault_sens_start, fault_sens_end,
                                                      fault_type, fault_functionpar, fault_a1, fault_a2)
                    level_pres = [round(elem, decimal_size) for elem in level_pres]
                    total_levels[node_id] = level_pres

            for j in range(0, self.wn.num_links):
                link_id = self.wn.link_name_list[j]
                if link_id not in flow_sensors:
                    continue
                flows = self.results.link['flowrate'][link_id]
                flows = flows[:len(self.time_stamp)]
                flows = [elem * 3600 for elem in flows]
                if f'flow{link_id}' in fault_objid.keys():
                    index_id = fault_index_ids[f'flow{link_id}']
                    flows = self.sensorfaultmodels(flows, index_id, fault_sens_start, fault_sens_end,
                                                      fault_type, fault_functionpar, fault_a1, fault_a2)
                flows = [round(elem, decimal_size) for elem in flows]
                total_flows[link_id] = flows

            # Create a Pandas dataframe from the data.
            df1 = pd.DataFrame(total_pressures)
            df2 = pd.DataFrame(total_demands)
            df3 = pd.DataFrame(total_flows)
            df4 = pd.DataFrame(total_levels)
            # Create a Pandas Excel writer using XlsxWriter as the engine.
            writer = pd.ExcelWriter(f'{results_folder}scenario{str(scenario)}\\Measurements.xlsx', engine='xlsxwriter')

            # Convert the dataframe to an XlsxWriter Excel object.
            # Pressures (m), Demands (m^3/h), Flows (m^3/h), Levels (m)
            df1.to_excel(writer, sheet_name='Pressures (m)', index=False)
            df2.to_excel(writer, sheet_name='Demands (L_h)', index=False)
            df3.to_excel(writer, sheet_name='Flows (m3_h)', index=False)
            df4.to_excel(writer, sheet_name='Levels (m)', index=False)

            # Close the Pandas Excel writer and output the Excel file.
            writer.save()

            try:
                os.remove('self.wn.pickle_fault')
            except:
                pass
        else:
            print('Results empty.')
            logging.info('Results empty.')
            return -1


if __name__ == '__main__':

    # Create tic / toc
    t = time.time()

    # Call dataset creator
    L = DatasetCreator()
    # Create scenario one-by-one
    scenarios = list(set(sensor_faults['scenario']))

    for scenario in scenarios:
        L.dataset_generator(scenario)

    print(f"Dataset completed.")
    logging.info(f"Dataset completed.")
    print(f'Total Elapsed time is {str(time.time() - t)} seconds.')
    logging.info(f'Total Elapsed time is {str(time.time() - t)} seconds.')
    os.startfile(logfilename)