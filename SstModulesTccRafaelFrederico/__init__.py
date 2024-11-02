import threading
import time as t
from termcolor import colored
from sst_reader import run_reader
from sst_feeder import run_feeder
from sst_data_capture import run_data_capture


RBD = {"SST_FILE_LOCATION": "", "SST_FILE_NAME": "", "SST_FILE_TYPE": "", "SST_TYPE": "", "SST_PREFIX": "", "SST_DATE": "",
"SST_TIME": "", "SST_PYTHON_DATA": "", "SST_VARS_NAMES": [], "SST_VARS_SIZES": [], "SST_VARS_RANGES": [],
"SST_PYTHON_BINARY_FORMAT": [], "SST_VARS_TOTAL_SIZE": 0}

XML_INFOS = {"XML_LOCATION": "", "XML_FILE_NAME": ""}

READ_MEMORY_DATA = {}

ring_list = list()


def get_file_type(file, file_location):
    RBD["SST_FILE_NAME"] = file
    RBD["SST_FILE_LOCATION"] = file_location + file
    
    if("rs" in file or "rf" in file):
        file_datime = file.split(".")
        prefix_date = file_datime[0]
        suffix_time = file_datime[1]
    else:
        prefix_date = file
        suffix_time = ""
        
    
    if("rs" in file):
        RBD["SST_TYPE"] = "Integration"
        RBD["SST_FILE_TYPE"] = "Data"
        RBD["SST_PREFIX"] = "rs"
    elif("rf" in file):
        RBD["SST_TYPE"] = "Subintegration"
        RBD["SST_FILE_TYPE"] = "Data"
        RBD["SST_PREFIX"] = "rf"
    elif("bi" in file):
        RBD["SST_TYPE"] = "Auxiliary"
        RBD["SST_FILE_TYPE"] = "Auxiliary"
        RBD["SST_PREFIX"] = "bi"
        


def main():
    file_name = "rs1170906.1600"
    file_location = "DATA/"
    xml_name = "SSTDataFormatTimeSpanTable.xml"
    xml_file_location = "XML_TABLES/"
    
    get_file_type(file_name, file_location)

    if(RBD["SST_TYPE"] == "Integration"):
        clock = 40
    elif(RBD["SST_TYPE"] == "Subintegration"):
        clock = 5
    elif(RBD["SST_TYPE"] == "Auxiliary"):
        clock = 1000

    print(colored("\nStarting reader module...", "green",  attrs = ["bold"]))
    thread_reader = threading.Thread(target=run_reader, args=(RBD, XML_INFOS, READ_MEMORY_DATA, file_name, xml_name, file_location, xml_file_location, clock))
    thread_reader.start()

    print(colored(80*"-", "red",  attrs = ["bold"]))
    t.sleep((clock*5)/1000)

    print(colored("\nStarting feeder module...", "green",  attrs = ["bold"]))
    thread_feeder = threading.Thread(target=run_feeder, args=(READ_MEMORY_DATA, RBD["SST_TYPE"], ring_list, clock))
    thread_feeder.start()

    print(colored(80*"-", "red",  attrs = ["bold"]))
    t.sleep((clock*5)/1000)

    print(colored("\nStarting data_capture module...", "green",  attrs = ["bold"]))
    thread_data_capture = threading.Thread(target=run_data_capture, args=(RBD["SST_TYPE"], RBD["SST_DATE"], RBD["SST_TIME"], ring_list, clock))
    thread_data_capture.start()


if __name__ == "__main__":
    main()
