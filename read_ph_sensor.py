import serial
import time
import sched
import csv
from datetime import datetime
import os
import re

def read_instrument(port, baudrate, timeout=2):
    with serial.Serial(port, baudrate, timeout=timeout) as ser:
        while True:
            ser.write(b'\r')
            print("wakeup!")
            time.sleep(0.5)  # Short delay between attempts
            response = ser.readline().decode('ascii').strip()
            print(response)
            if response.startswith("NAK"):
                break
        
        # Send the TS command
        ser.write(b'ts\r')
        print("Sent TS command")
        
        # Read lines until we get one starting with '#'
        while True:
            response = ser.readline().decode('ascii').strip()
            print(response)
            if response.startswith('#'):
                return response

def parse_data(data):
    # Use regex to split the data, handling variable whitespace
    pattern = r'#(\d+)\s+(\d{2}/\d{2}/\d{4})\s+(\d{2}:\d{2}:\d{2})\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)'
    match = re.match(pattern, data)
    
    if match:
        samp_num = int(match.group(1))
        datetime_str = f"{match.group(2)} {match.group(3)}"
        values = [float(match.group(i)) for i in range(4, 15)]
        
        return [samp_num, datetime_str] + values
    else:
        raise ValueError("Invalid data format")

def log_data(filename, data):
    file_exists = os.path.isfile(filename)
    
    with open(filename, 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        
        if not file_exists:
            # Write header if file doesn't exist
            csvwriter.writerow(['samp_num', 'datetime', 'v_bat', 'v_bias_pos', 'v_bias_neg', 
                                't_board', 'h_board', 'vrse', 'vrse_std', 'cevk', 'cevk_std', 
                                'ce_ik', 'i_sub'])
        
        csvwriter.writerow(data)

def scheduled_reading(scheduler, port, baudrate, filename):
    try:
        raw_data = read_instrument(port, baudrate)
        
        if raw_data:
            parsed_data = parse_data(raw_data)
            log_data(filename, parsed_data)
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
    scheduler.enter(10, 1, scheduled_reading, (scheduler, port, baudrate, filename))

if __name__ == "__main__":
    PORT = '/dev/tty.usbserial-FT9439MT0'  # Adjust this to your serial port
    BAUDRATE = 115200  
    LOGFILE = 'instrument_data.csv'  # Name of the CSV file
    
    s = sched.scheduler(time.time, time.sleep)
    
    # schedule first reading immediately
    s.enter(0, 1, scheduled_reading, (s, PORT, BAUDRATE, LOGFILE))
    
    print(f"Starting scheduled readings every 10 seconds. Logging to {LOGFILE}. Press Ctrl+C to stop.")
    
    try:
        s.run()
    except KeyboardInterrupt:
        print("\nScheduled readings stopped.")