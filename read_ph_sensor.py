import serial
import time
import sched
import csv
from datetime import datetime
import os
import re
import sqlite3
from calc_pH_DeepSeapHOx import calc_pH

READ_TIME = 10 # time between readings
PORT = 'COM6'  # Adjust this to your serial port
#PORT = '/dev/tty.usbserial-FT9439MT0'  # Adjust this to your serial port
BAUDRATE = 115200  
# LOGFILE = 'C:/Users/CSL 2/Documents/LOCNESS_data/pH_data.csv'  # Name of the CSV file
# DB_PATH = 'C:/Users/CSL 2/Documents/LOCNESS_data/data.db' # Path to SQLite DB
LOGFILE = 'pH_data.csv'  # Name of the CSV file
DB_PATH = 'data.db' # Path to SQLite DB

# Fixed values for rough pH calibration
TEMP = 25
SAL = 30
K0 = -1.42765
K2 = -0.001037
    
def read_instrument(port, baudrate, timeout=2):
    with serial.Serial(port, baudrate, timeout=timeout) as ser:
        print("wakeup!")
        while True:
            ser.write(b'\r')
            bytesToRead = ser.in_waiting
            response = ser.read(bytesToRead).decode('ascii').strip()
            if 'NAK' in response:
                print(response)
                break
            time.sleep(0.3)  # Short delay between attempts
        
        # Send the TS command
        ser.write(b'ts\r')
        print("Sent TS command")
        
        # Read lines until we get one starting with '#'
        # Need a t
        count = 0
        while True and count <= 20:
            response = ser.readline().decode('ascii').strip()
            #print(response)
            if response.startswith('#'):
                return response
            count+=1
            time.sleep(0.1)
        return 0
        
def parse_data(data, temp, sal, k0, k2):
    # Use regex to split the data, handling variable whitespace
    pattern = r'#(\d+)\s+(\d{2}/\d{2}/\d{4})\s+(\d{2}:\d{2}:\d{2})\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)'
    match = re.match(pattern, data)
    
    if match:
        samp_num = int(match.group(1))
        datetime_str = f"{match.group(2)} {match.group(3)}"
        values = [float(match.group(i)) for i in range(4, 15)]
        timestamp = time.time()
        ts = datetime.fromtimestamp(timestamp)
        ph_free, ph_tot = calc_pH(values[8], 0, temp, sal, k0, k2)
        return [ts, samp_num, datetime_str] + values + [temp, sal, k0, k2, ph_free, ph_tot]
    else:
        raise ValueError("Invalid data format")


def log_data_db(data):
    var_string = ', '.join('?' * len(data))
    query_string = f"INSERT INTO ph VALUES ({var_string});"
    c.execute(query_string, data)
    conn.commit()

def log_data(filename, data):
    file_exists = os.path.isfile(filename)
    
    with open(filename, 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        
        if not file_exists:
            # Write header if file doesn't exist
            csvwriter.writerow(['pc_time', 'samp_num', 'ph_time', 'v_bat', 'v_bias_pos', 'v_bias_neg', 
                                't_board', 'h_board', 'vrse', 'vrse_std', 'cevk', 'cevk_std', 
                                'ce_ik', 'i_sub', 'cal_temp', 'cal_sal', 'k0', 'k2', 'pH_free', 'ph_total'])
        
        csvwriter.writerow(data)

def scheduled_reading(scheduler, port, baudrate, filename):
    try:
        raw_data = read_instrument(port, baudrate)
        
        if raw_data:
            parsed_data = parse_data(raw_data, TEMP, SAL, K0, K2)
            log_data(filename, parsed_data)
            log_data_db(parsed_data)
            print("Logged data:", parsed_data)
        else:
            print(f"Time: {datetime.now().strftime('%m/%d/%Y %H:%M:%S')}")
            print("No data received from the instrument")
            print()
    
    except serial.SerialException as e:
        print(f"Error: {e}")
        print("Please check your serial port configuration and connections.")
        print()
    
    # Schedule the next reading
    scheduler.enter(READ_TIME, 1, scheduled_reading, (scheduler, port, baudrate, filename))

if __name__ == "__main__":

    # Connect to the SQLite database 
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()

    # Create a table to store the data
    # or append if it already exists
    c.execute('''CREATE TABLE IF NOT EXISTS ph
            (samp_num INTEGER, ph_timestamp TEXT, 
             v_bat REAL, v_bias_pos REAL, v_bias_neg REAL, 
             t_board REAL, h_board REAL, vrse REAL, vrse_std REAL, 
             cevk REAL, cevk_std REAL, ce_ik REAL, i_sub REAL,
             cal_temp REAL, cal_sal REAL, k0 REAL, k2 REAL,
             pH_free REAL, ph_total REAL
             )''')

    s = sched.scheduler(time.time, time.sleep)
    
    # schedule first reading immediately
    s.enter(0, READ_TIME, scheduled_reading, (s, PORT, BAUDRATE, LOGFILE))
    
    print(f"Starting scheduled readings every 10 seconds. Logging to {LOGFILE}. Press Ctrl+C to stop.")
    
    try:
        s.run()
    except KeyboardInterrupt:
        conn.close()
        print("\nScheduled readings stopped.")