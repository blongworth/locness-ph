import serial
import time
import sched
import csv
from datetime import datetime
import os
import re
import sqlite3
import yaml
from calc_pH_DeepSeapHOx import calc_pH
import logging
#from datetime import timezone

# Read the configuration file
with open('config.yaml', 'r') as file:
    config = yaml.safe_load(file)

# Access configuration values
READ_TIME = config['read_time']
PORT = config['sensor']['port']
BAUDRATE = config['sensor']['baudrate']
DATAFILE = config['file']['data']
LOGFILE = config['file']['log']
DB_PATH = config['file']['db']
TEMP = config['calibration']['temp']
SAL = config['calibration']['sal']
K0 = config['calibration']['k0']
K2 = config['calibration']['k2']

# Set the timezone to UTC
#logging.Formatter.converter = lambda *args: datetime.now(timezone.utc).timetuple()
# Configure logging
logging.basicConfig(
    # filename=LOGFILE,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%dT%H:%M:%S%z',
     handlers=[
        logging.FileHandler(LOGFILE),
        logging.StreamHandler()
    ]
)

# Start logger
logger = logging.getLogger(__name__)

def read_instrument(port, baudrate, timeout=2):
    with serial.Serial(port, baudrate, timeout=timeout) as ser:
        logger.info("Wake mFET")
        while True:
            ser.write(b"\r")
            bytesToRead = ser.in_waiting
            response = ser.read(bytesToRead).decode("ascii").strip()
            if "NAK" in response:
                logger.info(f"Wake response: {response}")
                break
            time.sleep(0.3)  # Short delay between attempts

        # Send the TS command
        ser.write(b"ts\r")
        logger.info("Sent TS command")

        # Read lines until we get one starting with '#'
        count = 0
        while True and count <= 20:
            response = ser.readline().decode("ascii").strip()
            # print(response)
            if response.startswith("#"):
                return response
            count += 1
            time.sleep(0.1)
        return 0

def parse_data(data, temp, sal, k0, k2):
    # Use regex to split the data, handling variable whitespace
    pattern = r"#(\d+)\s+(\d{2}/\d{2}/\d{4})\s+(\d{2}:\d{2}:\d{2})\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    match = re.match(pattern, data)

    if match:
        samp_num = int(match.group(1))
        datetime_str = f"{match.group(2)} {match.group(3)}"
        values = [float(match.group(i)) for i in range(4, 15)]
        timestamp = time.time()
        ts = datetime.fromtimestamp(timestamp)
        ph_free, ph_tot = calc_pH(values[5], 0, temp, sal, k0, k2, 0)
        return (
            [ts, samp_num, datetime_str] + values + [temp, sal, k0, k2, ph_free, ph_tot]
        )
    else:
        raise ValueError("Invalid data format")


def log_data_db(data):
    var_string = ", ".join("?" * len(data))
    query_string = f"INSERT INTO ph VALUES ({var_string});"
    c.execute(query_string, data)
    conn.commit()


def log_data(filename, data):
    file_exists = os.path.isfile(filename)

    with open(filename, "a", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)

        if not file_exists:
            # Write header if file doesn't exist
            csvwriter.writerow(
                [
                    "pc_time",
                    "samp_num",
                    "ph_time",
                    "v_bat",
                    "v_bias_pos",
                    "v_bias_neg",
                    "t_board",
                    "h_board",
                    "vrse",
                    "vrse_std",
                    "cevk",
                    "cevk_std",
                    "ce_ik",
                    "i_sub",
                    "cal_temp",
                    "cal_sal",
                    "k0",
                    "k2",
                    "ph_free",
                    "ph_total",
                ]
            )

        csvwriter.writerow(data)


def scheduled_reading(scheduler, port, baudrate, filename):
    try:
        raw_data = read_instrument(port, baudrate)

        if raw_data:
            parsed_data = parse_data(raw_data, TEMP, SAL, K0, K2)
            log_data(filename, parsed_data)
            log_data_db(parsed_data)
            logger.info(f"Logged data: {parsed_data}")
        else:
            logger.error("No data received from the instrument")

    except serial.SerialException as e:
        logger.error(e)
        exit(1)

    # Schedule the next reading
    scheduler.enter(
        READ_TIME, 1, scheduled_reading, (scheduler, port, baudrate, filename)
    )


if __name__ == "__main__":
    # Connect to the SQLite database
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()

    # Create a table to store the data
    # or append if it already exists
    c.execute("""CREATE TABLE IF NOT EXISTS ph
            (pc_timestamp TEXT, samp_num INTEGER, ph_timestamp TEXT, 
             v_bat REAL, v_bias_pos REAL, v_bias_neg REAL, 
             t_board REAL, h_board REAL, vrse REAL, vrse_std REAL, 
             cevk REAL, cevk_std REAL, ce_ik REAL, i_sub REAL,
             cal_temp REAL, cal_sal REAL, k0 REAL, k2 REAL,
             ph_free REAL, ph_total REAL
             )""")

    s = sched.scheduler(time.time, time.sleep)

    # schedule first reading immediately
    s.enter(0, READ_TIME, scheduled_reading, (s, PORT, BAUDRATE, DATAFILE))

    logger.info(f"Starting scheduled readings every {READ_TIME} seconds. Logging to {DATAFILE}. Press Ctrl+C to stop.")

    try:
        s.run()
    except KeyboardInterrupt:
        conn.close()
        logger.info("Scheduled readings stopped.")
