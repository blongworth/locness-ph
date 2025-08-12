import csv
import logging
import os
import sys
import re
import sched
import serial
import sqlite3
import time
from datetime import datetime
import yaml
from isfetphcalc import calc_ph

# TODO: use integer timestamp instead of datetime for performance

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
    datefmt='%Y-%m-%dT%H:%M:%S',
     handlers=[
        logging.FileHandler(LOGFILE),
        logging.StreamHandler()
    ]
)

# Start logger
logger = logging.getLogger(__name__)

def read_instrument(port, baudrate, timeout=2, polled=True):
    with serial.Serial(port, baudrate, timeout=timeout) as ser:
        if polled == True:
            logger.debug("Wake mFET")
            while True:
                ser.write(b"\r")
                bytesToRead = ser.in_waiting
                try:
                    response = ser.readline(bytesToRead).decode("ascii").strip()
                except UnicodeDecodeError as e:
                    logger.error(e)
                if "NAK" in response:
                    logger.debug(f"Wake response: {response}")
                    time.sleep(0.1)  # Short delay between attempts
                    break
                time.sleep(0.1)  # Short delay between attempts

            # Send the TS command
            ser.write(b"ts\r")
            logger.debug("Sent TS command")

        # Read lines until we get one starting with '#'
        count = 0
        while True and count <= 20:
            try:
                response = ser.readline().decode("ascii").strip()
            except UnicodeDecodeError as e:
                logger.error(e)
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
        ph_free, ph_tot = calc_ph(values[5], 0, temp, sal, k0, k2, 0)
        ph_free = round(float(ph_free), 4)
        ph_tot = round(float(ph_tot), 4)
        return {
            "datetime_utc": ts,
            "samp_num": samp_num,
            "ph_timestamp": datetime_str,
            "v_bat": values[0],
            "v_bias_pos": values[1],
            "v_bias_neg": values[2],
            "t_board": values[3],
            "h_board": values[4],
            "vrse": values[5],
            "vrse_std": values[6],
            "cevk": values[7],
            "cevk_std": values[8],
            "ce_ik": values[9],
            "i_sub": values[10],
            "cal_temp": temp,
            "cal_sal": sal,
            "k0": k0,
            "k2": k2,
            "ph_free": ph_free,
            "ph_total": ph_tot
        }
    else:
        logger.error(f"Data format error: {data}")
        raise ValueError("Invalid data format")

def log_data_db(data):
    # convert the first element of data from datetime to integer timestamp
    data = dict(data)
    if isinstance(data["datetime_utc"], datetime):
        data["datetime_utc"] = int(data["datetime_utc"].timestamp())
    columns = [
        "datetime_utc",
        "samp_num",
        "ph_timestamp",
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
        "ph_total"
    ]
    values = [data[col] for col in columns]
    col_string = ", ".join(columns)
    var_string = ", ".join(["?" for _ in columns])
    query_string = f"INSERT INTO ph ({col_string}) VALUES ({var_string});"
    try:
        c.execute(query_string, values)
        conn.commit()
    except Exception as e:
        logger.error(f"Database error: {e}")
        conn.rollback()

def log_data(filename, data):
    file_exists = os.path.isfile(filename)
    data = dict(data)
    values = [
        data["datetime_utc"],
        data["samp_num"],
        data["ph_timestamp"],
        data["v_bat"],
        data["v_bias_pos"],
        data["v_bias_neg"],
        data["t_board"],
        data["h_board"],
        data["vrse"],
        data["vrse_std"],
        data["cevk"],
        data["cevk_std"],
        data["ce_ik"],
        data["i_sub"],
        data["cal_temp"],
        data["cal_sal"],
        data["k0"],
        data["k2"],
        data["ph_free"],
        data["ph_total"]
    ]
    with open(filename, "a", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        if not file_exists:
            # Write header if file doesn't exist
            csvwriter.writerow([
                "datetime_utc",
                "samp_num",
                "ph_timestamp",
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
            ])
        csvwriter.writerow(values)

def scheduled_reading(scheduler, port, baudrate, filename):
    # Schedule the next reading
    scheduler.enter(
        READ_TIME, 1, scheduled_reading, (scheduler, port, baudrate, filename)
    )

    # If read time is 0, sensor is in internally timed mode.
    # Don't poll sensor, and wait for next reading as soon as
    # previous is completed
    if READ_TIME == 0:
        polled = False
    else:
        polled = True

    try:
        raw_data = read_instrument(port, baudrate, polled=polled)

        if raw_data:
            parsed_data = parse_data(raw_data, TEMP, SAL, K0, K2)
            log_data(filename, parsed_data)
            log_data_db(parsed_data)
            # Convert datetime (first element) to ISO-8601 string for logging
            logged_data = dict(parsed_data)
            if isinstance(logged_data["datetime_utc"], datetime):
                logged_data["datetime_utc"] = logged_data["datetime_utc"].replace(microsecond=0).isoformat()
            logger.debug(f"Logged data: {logged_data}")
            logger.info(
                f"vm: {parsed_data['v_bat']:.3g}, vb: {parsed_data['v_bias_pos']:.3g}, vrse: {parsed_data['vrse']:.3g}, vrsd: {parsed_data['vrse_std']:.3g}, ph={parsed_data['ph_total']:.3g}"
            )
        else:
            logger.error("No data received from the instrument")

    except serial.SerialException as e:
        logger.error(e)
        exit(1)

def ensure_database_ready(db_path):
    """Quick check that database is properly initialized"""
    try:
        conn = sqlite3.connect(DB_PATH)
        conn.execute('SELECT 1 FROM ph LIMIT 1')
        conn.close()
        return True
    except sqlite3.OperationalError:
        return False

def main():
    """Main function to run the scheduled readings"""
    if not ensure_database_ready(DB_PATH):
        logger.error("Database not initialized. Set up with locness-datamanager first.")
        sys.exit(1)
    
    # Connect to the SQLite database
    global conn, c
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()

    s = sched.scheduler(time.time, time.sleep)

    # schedule first reading immediately
    s.enter(0, READ_TIME, scheduled_reading, (s, PORT, BAUDRATE, DATAFILE))

    logger.info(f"Starting scheduled readings every {READ_TIME} seconds. Logging to {DATAFILE}. Press Ctrl+C to stop.")

    try:
        s.run()
    except KeyboardInterrupt:
        conn.close()
        logger.info("Scheduled readings stopped.")

if __name__ == "__main__":
    main()