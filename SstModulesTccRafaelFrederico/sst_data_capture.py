import time as t
import traceback
from lxml import etree
from datetime import datetime
from termcolor import colored
from email.generator import Generator
from email.mime.text import MIMEText as MT


def extract_ring_list_values(sst_type, sst_date, sst_time, ring_list, clock):
    index = 0
    count_errors = 0
    ring_size = ring_list[index][7]
    root = etree.Element("SSTDataSet", {"DataType":sst_type})

    while True:
        try:
            start_time = int(round(t.time() * 1000))

            if(("end" in ring_list) and (ring_list.index("end") == index) or (count_errors == 10)):
                print(colored("Stopping data capture...\n", "green", attrs = ["bold"]))
                break

            if(index == ring_size):
                index = 0

            if(ring_list[index]):
                sdmDataSubset = etree.Element("SdmDataSubset",{
                                                "time": str(ring_list[index][0]),
                                                "target": str(ring_list[index][8]),
                                                "opmode": str(ring_list[index][9])})

                adc0 = etree.Element("Adcval0")
                adc1 = etree.Element("Adcval1")
                adc2 = etree.Element("Adcval2")
                adc3 = etree.Element("Adcval3")
                adc4 = etree.Element("Adcval4")
                adc5 = etree.Element("Adcval5")

                adc0.text = str(ring_list[index][1])
                adc1.text = str(ring_list[index][2])
                adc2.text = str(ring_list[index][3])
                adc3.text = str(ring_list[index][4])
                adc4.text = str(ring_list[index][5])
                adc5.text = str(ring_list[index][6])

                sdmDataSubset.append(adc0)
                sdmDataSubset.append(adc1)
                sdmDataSubset.append(adc2)
                sdmDataSubset.append(adc3)
                sdmDataSubset.append(adc4)
                sdmDataSubset.append(adc5)

                root.append(sdmDataSubset)

                print(colored(f"Data Capture: extracted line: {index} \n values: {ring_list[index]} ", "green", attrs = ["bold"]))
                ring_list[index] = None
                count_errors = 0
                index += 1

            else:
                print(colored(f"No more data at line: {index} ", "green",  attrs = ["bold"]))
                count_errors += 1


            end_time = int(round(t.time() * 1000))
            sleep_time = ((clock - (end_time - start_time)) / 1000)
            t.sleep(sleep_time)

        except Exception:
            print(colored(f"Could not read the line: {index} from feeder, trying again...", "red",  attrs = ["bold"]))
            count_errors += 1
            traceback.print_exc()

    xmlstr = etree.tostring(root, encoding = "iso-8859-1")
    return xmlstr



def create_mime_file(xmlstr, sst_type, sst_date, sst_time):
    now_date = datetime.now()
    mime_body = MT(str(xmlstr))
    mime_body["mime_generate_date"] = str(now_date)
    mime_body["sst_observation_type"] = str(sst_type)
    mime_body["sst_observation_date"] = str(sst_date) + ":" + str(sst_time)
    filename = str(sst_type) + "_" + str(now_date) + ".txt"

    with open(filename, "w") as f:
        gen = Generator(f)
        gen.flatten(mime_body)



def run_data_capture(sst_type, sst_date, sst_time, ring_list, clock):
    xmlstr = extract_ring_list_values(sst_type, sst_date, sst_time, ring_list, clock)
    print(colored("Generating MIME table...\n", "green", attrs = ["bold"]))
    create_mime_file(xmlstr, sst_type, sst_date, sst_time)

