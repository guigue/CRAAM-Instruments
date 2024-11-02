import os
import sys
import json
import struct
import time as t
import numpy as np
from termcolor import colored
import xml.etree.ElementTree as xmltree
from datetime import datetime, timedelta


def get_file_datetime(RBD, file, file_location):
    RBD["SST_FILE_NAME"] = file
    RBD["SST_FILE_LOCATION"] = file_location + file
    
    if("rs" in file or "rf" in file):
        file_datime = file.split(".")
        prefix_date = file_datime[0]
        suffix_time = file_datime[1]
    else:
        prefix_date = file
        suffix_time = ""
    
    
    if(len(prefix_date) == 8):
        date = datetime(year=int(prefix_date[2:4]) + 1900, month=int(prefix_date[4:6]), day=int(prefix_date[6:8]))
        date_str = str(date).split(" ")[0]
        RBD["SST_DATE"] = date_str
        
    else:
        date = datetime(year=int(prefix_date[2:5]) + 1900, month=int(prefix_date[5:7]), day=int(prefix_date[7:9]))
        date_str = str(date).split(" ")[0]
        RBD["SST_DATE"] = date_str
    
    
    if(len(suffix_time) == 4):
        time = suffix_time[0:2] + ':' + suffix_time[2:4]
        RBD["SST_TIME"] = time
    elif(len(suffix_time) == 3):
        time = suffix_time[0:1] + ':' + suffix_time[1:3]
        RBD["SST_TIME"] = time
    else:
        time = "00:00"
        RBD["SST_TIME"] = time
        


def find_xml_description_file(RBD, XML_INFOS, xml_file, xml_file_location):
    xml_path = xmltree.parse(xml_file_location + xml_file)
    xml_tables = xml_path.getroot()
    
    for table in xml_tables:
        sst_file_type = table[0].text
        initial_date = table[1].text
        final_date = table[2].text
        xml_name = table[3].text
        
        if(table and RBD["SST_FILE_TYPE"] == sst_file_type):
            if(RBD["SST_DATE"] >= initial_date and RBD["SST_DATE"] <= final_date):
                XML_INFOS["XML_LOCATION"] = (xml_file_location + xml_name)
                XML_INFOS["XML_FILE_NAME"] = xml_name
                break

                

def get_xml_variables(RBD, XML_INFOS):
    xml_path = xmltree.parse(XML_INFOS["XML_LOCATION"])
    xml_tables = xml_path.getroot()
    
    python_data_type = get_sys_byte_order()
    vars_names = []
    vars_sizes = []
    vars_ranges = []
    vars_total_size = 0
    var_pointerdata = 0
    
    for table in xml_tables:
        var_name = table[0].text
        var_size = int(table[1].text)
        var_type = table[2].text
        
        vars_total_size += var_size
        
        vars_names.append(var_name)
        vars_sizes.append(var_size)
        vars_ranges.append([var_pointerdata, var_pointerdata + var_size - 1])
        
        var_pointerdata += var_size
        
        for i in range(var_size):
            
            if(var_type == "xs:int"):
                python_data_type += "i"
            elif(var_type == "xs:unsignedShort"): 
                python_data_type += "H"
            elif(var_type == "xs:short"):
                python_data_type += "h"
            elif(var_type == "xs:byte"):
                python_data_type += "b"
            elif(var_type == "xs:float"):
                python_data_type += "f"
                
        RBD["SST_VARS_NAMES"] = vars_names
        RBD["SST_VARS_SIZES"] = vars_sizes
        RBD["SST_VARS_RANGES"] = vars_ranges
        RBD["SST_PYTHON_DATA"] = python_data_type
        RBD["SST_VARS_TOTAL_SIZE"] = vars_total_size



def get_python_binary_data_format(RBD, XML_INFOS):
    xml_path = xmltree.parse(XML_INFOS["XML_LOCATION"])
    xml_tables = xml_path.getroot()
    
    python_binary_format = []
    
    for table in xml_tables:
        var_type = table[2].text
        
        if(var_type == "xs:int"):
            python_binary_format.append("int32")
        elif(var_type == "xs:unsignedShort"): 
            python_binary_format.append("uint16")
        elif(var_type == "xs:short"):
            python_binary_format.append("int16")
        elif(var_type == "xs:byte"):
            python_binary_format.append("byte")
        elif(var_type == "xs:float"):
            python_binary_format.append("float32")
        else:
            python_binary_format.append("float32")
    
    RBD["SST_PYTHON_BINARY_FORMAT"] = python_binary_format



def get_sys_byte_order():
    if(sys.byteorder == "little"):
        return "<"
    elif(sys.byteorder == "big"):
        return ">"
    else:
        return "="



def read_sst_raw_binary_data(RBD, XML_INFOS, READ_MEMORY_DATA, clock):
    python_data_type = RBD["SST_PYTHON_DATA"]
    python_binary_format = RBD["SST_PYTHON_BINARY_FORMAT"]
    vars_names = RBD["SST_VARS_NAMES"]
    vars_ranges = RBD["SST_VARS_RANGES"]
    vars_sizes = RBD["SST_VARS_SIZES"]
    num_vars = len(vars_names)
    actual_read_line = 0
    
    if(os.path.isfile(RBD["SST_FILE_LOCATION"]) and os.access(RBD["SST_FILE_LOCATION"], os.R_OK)):
        file_description = os.open(RBD["SST_FILE_LOCATION"], os.O_RDONLY)
        file_mem_size = os.fstat(file_description).st_size
        table_mem_size = struct.calcsize(python_data_type)
        py_mem_size = int(file_mem_size / table_mem_size)
        
        for i in range(len(python_binary_format)):
            if(python_binary_format[i] == "int32"):
                nump_data = np.int32
            elif(python_binary_format[i] == "uint16"):
                nump_data = np.uint16
            elif(python_binary_format[i] == "int16"):
                nump_data = np.int16
            elif(python_binary_format[i] == "byte"):
                nump_data = np.byte
            elif(python_binary_format[i] == "float32"):
                nump_data = np.float32

            if(vars_sizes[i] == 1):
                mem_block = {vars_names[i]: np.array(np.empty([py_mem_size], nump_data))}
                READ_MEMORY_DATA.update(mem_block)
            else:
                mem_block = {vars_names[i]: np.array(np.empty([py_mem_size, vars_sizes[i]], nump_data))}
                READ_MEMORY_DATA.update(mem_block)
                
        for j in range(py_mem_size):
            try:
                start_time = int(round(t.time() * 1000))
                read_hexa = os.read(file_description, table_mem_size)
                file_data = struct.unpack(python_data_type, read_hexa)
                actual_read_line += 1
                
            except Exception as e:
                print("Line:", actual_read_line , e)
                break;
                         
            try:
                for z in range (num_vars):
                    actual_range = vars_ranges[z]
                    
                    if(actual_range[0] == actual_range[1]):
                        READ_MEMORY_DATA[vars_names[z]][j] = file_data[actual_range[0]]
                    else:
                        READ_MEMORY_DATA[vars_names[z]][j, ...] = file_data[actual_range[0] : actual_range[1]+1]
                        
                end_time = int(round(t.time() * 1000))
                sleep_time = ((clock - (end_time - start_time)) / 1000)
                t.sleep(sleep_time)
                
            except Exception as e:
                print("Unable to alocate variable", READ_MEMORY_DATA[vars_names[z]], "at file line", actual_read_line)
                error_flag = True
                print(e)
        
        percentage_rate = ((actual_read_line * 100) / py_mem_size)
        
        if(int(percentage_rate) == 100):
            show_json_data(RBD, XML_INFOS)
            print(colored(f"\nREADER: Lines successfully read: {actual_read_line} lines", "green", attrs = ["bold"]))
            print(colored("Stopping reader...\n", "green", attrs = ["bold"]))
        else:
            print("ERROR, possible to read", "%.2f" % percentage_rate, "% of the file")



def show_json_data(RBD, XML_INFOS):  
    RBD_JSON = json.dumps(RBD, indent = 2)
    XML_JSON = json.dumps(XML_INFOS, indent = 2)
    
    print(colored("\nRBD =", "blue",  attrs = ["bold"]), 
          colored(RBD_JSON, "blue"))
    
    print(colored("\nXML_INFOS =", "green", attrs = ["bold"]), 
          colored(XML_JSON, "green"))


            
def run_reader(RBD, XML_INFOS, READ_MEMORY_DATA, file_name, xml_name, file_location, xml_file_location, clock):
    get_file_datetime(RBD, file_name, file_location)
    find_xml_description_file(RBD, XML_INFOS, xml_name, xml_file_location)
    get_xml_variables(RBD, XML_INFOS)
    get_python_binary_data_format(RBD, XML_INFOS)
    read_sst_raw_binary_data(RBD, XML_INFOS, READ_MEMORY_DATA, clock)
    