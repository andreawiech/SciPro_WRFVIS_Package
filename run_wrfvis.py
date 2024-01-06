from wrfvis import cfg, core

# Specify the parameters
param = 'T'
lon = 11.3960  # Longitude for Innsbruck
lat = 47.2692  # Latitude for Innsbruck
zagl = 576.0   # Grid point elevation for Innsbruck
time_index = 3  # Time index for skewT plot
rad = 8000

# Specify the output directory
output_directory = cfg.output_directory

# Test the combined HTML generation function
combined_html = core.generate_combined_html(param, lon, lat, time_index, zagl, rad=rad, directory=output_directory)

# Print the path to the generated HTML file
print(f"Combined HTML file saved at: {combined_html}")