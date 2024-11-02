import time as t
import traceback
from termcolor import colored


def extract_dict_values(READ_MEMORY_DATA, sst_type, ring_list, clock):
    i = 0
    time_bin = 0
    ring_size = 0
    last_milli = 0
    count_errors = 0
    line_erros = []
    adcval_name = ""
    time = READ_MEMORY_DATA["time"]
    
    if(sst_type == "Auxiliary"):
        ring_size = 10
        time_bin = "1000 ms"
        adcval_name = "adc"
    elif(sst_type == "Integration"):
        ring_size = 100
        time_bin = "40 ms"
        adcval_name = "adcval"
    elif(sst_type == "Subintegration"):
        ring_size = 400
        time_bin = "5 ms"
        adcval_name = "adcval"
    
    while True:
        start_time = int(round(t.time() * 1000))
        if(count_errors >= 10):
            print(colored(f"\nFEEDER: Lines successfully extracted: {i} lines", "green",  attrs = ["bold"]))
            print(colored("Stopping feeder...\n", "green", attrs = ["bold"]))
            ring_list.append("end")
            check_data_erros(line_erros, i)
            break
        
        try: 
            actual_time = time[i]
            
            if(actual_time != 0 and last_milli != 0):
                limit_ms_error = int((time[i] - time[i - 1]) * 1.75)
                if((actual_time - last_milli) > limit_ms_error):
                    print(colored(f"Too long observation interval: \n"
                                  +f"expected: {time_bin} \n"
                                  +f"received: {(actual_time - last_milli)} ms",
                                  "red", attrs = ["bold"]))
            
            adcval0 = READ_MEMORY_DATA[adcval_name][i][0]
            adcval1 = READ_MEMORY_DATA[adcval_name][i][1]
            adcval2 = READ_MEMORY_DATA[adcval_name][i][2]
            adcval3 = READ_MEMORY_DATA[adcval_name][i][3]
            adcval4 = READ_MEMORY_DATA[adcval_name][i][4]
            adcval5 = READ_MEMORY_DATA[adcval_name][i][5]
            target = READ_MEMORY_DATA["target"][i]
            opmode = READ_MEMORY_DATA["opmode"][i]

            tuple = (actual_time, adcval0, adcval1, adcval2,
                    adcval3, adcval4, adcval5, ring_size,
                    target, opmode)
            
            save_to_ring_list(ring_list, ring_size, tuple, i)
            
            last_milli = actual_time
            i += 1
            
        except IndexError:
            line_erros.append(i)
            count_errors += 1

        end_time = int(round(t.time() * 1000))
        sleep_time = ((clock - (end_time - start_time)) / 1000)
        t.sleep(sleep_time)


            
def save_to_ring_list(ring_list, ring_size, tuple, i):
    if(i < ring_size):
        ring_list.append(tuple)
    else:
        ring_list[i%ring_size] = tuple
        

        
def check_data_erros(line_erros, last_line):
    size_erros = len(line_erros)
    
    if(line_erros[0] and line_erros[0] != last_line):
        print(colored("Could not read the lines: ", "red",  attrs = ["bold"]))
        
        for i in range (0, size_erros):
            print(line_erros[i])
    else:
        pass
        
        
    
def run_feeder(READ_MEMORY_DATA, sst_type, ring_list, clock):
    try:
        extract_dict_values(READ_MEMORY_DATA, sst_type, ring_list, clock)
    except Exception:
        traceback.print_exc()
        