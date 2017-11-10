import visa
import csv
import numpy
from time import sleep

global sweep_data

def afg_set_freq(freq_hz):
    afg.write("SOURce1:FREQuency:CW " + str(freq_hz) + " Hz")

def afg_set_amp(amp_v):
    afg.write("SOURce1:VOLTage:LEVel:IMMediate:HIGH " + str(amp_v / 2.0) + "V")
    afg.write("SOURce1:VOLTage:LEVel:IMMediate:LOW " + str(amp_v / -2.0) + "V")
    afg.write("SOURce1:VOLTage:LEVel:IMMediate:OFFSet " + str(0) + "V")

def import_sweep_table(file):
    global sweep_data
    with open(file, 'r') as f:
        reader = csv.reader(f)
        next(reader, None) # Skip headers
        sweep_data = list(reader)
    f.close()

def run_sweep():
    global sweep_data
    for x in sweep_data:
        x = map(float, x)
        freq_list = numpy.arange(x[1], x[2]+x[3], x[3])
        afg_set_amp(x[0])  # Set Amplitude
        print "Setting Amplitude to " + str(x[0]) + "V"
        print "Dwell time is " + str(x[4]) + "s"
        for freq in freq_list:
            print "Setting Frequency to " + str(freq) + "Hz"
            afg_set_freq(freq) # Set Frequency
            sleep(x[4])

def calc_sweep_time():
    global sweep_data
    total_time = 0
    for x in sweep_data:
        x = map(float, x)
        freq_list = numpy.arange(x[1], x[2]+x[3], x[3])
        for freq in freq_list:
            total_time += x[4]
    m, s = divmod(total_time, 60)
    h, m = divmod(m, 60)
    print "Total Test Time: %d:%02d:%02d" % (h, m, s)

if __name__ == "__main__":
    rm = visa.ResourceManager()

    global sweep_data
    import_sweep_table('sweep_table.csv')

    calc_sweep_time()
    raw_input("Press enter to continue.")

    try:
        afg = rm.open_resource('USB0::0x0699::0x0352::C010959::INSTR')
    except:
        print "Couldn't locate instrument. Available instruments:"
        print(rm.list_resources())
        exit()

    run_sweep()
    afg_set_amp(0)  # Set Amplitude
    exit()