# Default configuration
# Local configuration should be in config.py

READ_TIME = 10  # time between readings
PORT = "COM6"  # Adjust this to your serial port
# PORT = '/dev/tty.usbserial-FT9439MT0'  # Adjust this to your serial port
BAUDRATE = 115200
# LOGFILE = 'C:/Users/CSL 2/Documents/LOCNESS_data/pH_data.csv'  # Name of the CSV file
# DB_PATH = 'C:/Users/CSL 2/Documents/LOCNESS_data/data.db' # Path to SQLite DB
LOGFILE = "pH_data.csv"  # Name of the CSV file
DB_PATH = "data.db"  # Path to SQLite DB

# Fixed values for rough pH calibration
TEMP = 25
SAL = 30
K0 = -1.42765
K2 = -0.001037