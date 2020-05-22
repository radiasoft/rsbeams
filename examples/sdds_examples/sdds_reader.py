from rsbeams.rsdata.SDDS import readSDDS

sdds_file = './support/example_sdds_file.sig'

# Initialize the reader
reader = readSDDS(sdds_file)
# Read file contents
reader.read()

# Data of different formats (columns and parameters) may be accessed by attributes of the same name
columns = reader.columns
parameters = reader.parameters

# Data is stored in structured NumPy arrays
# Array shape is based on the number of data entries and pages
# field shape: (page, entries)
print("Each column field has shape: {} pages and {} entries".format(*columns.shape))
# parameters can only store one entry each
print("Each parameter field has shape: {} pages and {} entries\n".format(*parameters.shape))


# To access each field you can use the usual structured NumPy array syntax
print("Data from the 'Sx' field is accessed by `columns['Sx']`: {}\n".format(columns['Sx'][0, :10]))
print("You can get a list of all fields with the `dtype` attribute.\nparameters.dtype={} ".format(parameters.dtype))