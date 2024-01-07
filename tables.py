"""
Module for generating HTML tables with skiing condition assessments based on WRF model parameters.

Author
------
Joep van Noort, Malte Hildebrandt

Functions
---------
1. `get_match_value(match)`: 
    Extracts output variables from `snowcheck.snow_variables` for each integer match in 'rows' HTML string. 
    Used in `weather_table`.

2. `weather_table(lon, lat, ds)`: 
    Generates an HTML file with a table indicating skiing conditions. 
    Uses `snowcheck.snow_variables` for assessing WRF parameter values with dynamic coloring.

3. `geographical_table(lon, lat, ds)`: 
    Generates an HTML file with a table providing location-dependent information for potential winter skiing. 
    Uses `snowcheck.mountain_check` with dynamic coloring.

4. `html_page(html_table, html_table2, filename="tables.html")`: 
    Combines two HTML tables into a file saved to the specified directory.

Dependencies
------------
- numpy as np
- pandas as pd
- xarray as xr
- re
- wrfvis.cfg
- snowcheck

Usage
-----
1. Import the module: `import wrfvis.tables.
2. Call the desired functions based on skiing condition assessments.
3. Ensure that the required dependencies are installed before using the module.
"""
import os
import pandas as pd
import xarray as xr
import re

from wrfvis import cfg, snowcheck


def get_match_value(match):
    """
    @author: Joep van Noort

    This function is called within the weather_table function as a regex
    argument. It works similar to a loop, where a counter is used to
    sequentially obtain output variables from the snowcheck.snow_variables
    function for each integer match in the 'rows' html-string. This allows
    the dynamically generated colored html table within the snow_variables
    function to retain its colors (based on a numerical value) and change
    the displayed values in the final resulting html file to natural language.

    Parameters
    ----------
    match : match to predefined regex pattern

    Returns
    -------
    output : natural language or numerical value corresponding to the
             snow_variables function output indeces 1, 3, 5, 7.

    """
    global counter
    global timestep
    output = snowcheck.snow_variables(lon, lat, ds, timestep)[counter]
    counter += 2
    if counter > 7:
        counter = 1
        timestep += 1
    if timestep > 35:
        timestep = 0
    return output


def weather_table(lon, lat, ds):
    """
    @authors: Malte Hildebrandt, Joep van Noort

    This function generates an html file containing a table in which certain
    (modified) WRF parameter values provide information about weather
    conditions for skiing.

    The function that is used to obtain this information is the snow_variables
    function from the snowcheck package. The table is then dynamically
    generated with value-dependent coloring.

    Parameters
    ----------
    lon : user input longitude value in degrees East
    lat : user input latitude value in degrees North
    ds : wrfout dataset
    time : WRF model timestep (0-35, default 24)

    Returns
    -------
    html_table : html table in string format. A html file is generated in
                 the backgorund.

    """
    # create timeseries dataset for weather variables
    timestep = range(36)
    time_list = []
    snowcover_list = []
    snowfall_list = []
    sun_list = []
    wind_list = []

    for i in timestep:
        date = ds['XTIME'][i].values
        time_value = date.astype(str)[5:10] + ' ' + date.astype(str)[11:16]
        time_list.append(time_value)
        snowcover_list.append(snowcheck.snow_variables(lon, lat, ds, i)[0])
        snowfall_list.append(snowcheck.snow_variables(lon, lat, ds, i)[2])
        sun_list.append(snowcheck.snow_variables(lon, lat, ds, i)[4])
        wind_list.append(snowcheck.snow_variables(lon, lat, ds, i)[6])

    data = {'Time': time_list,
            'Snowcover': snowcover_list,
            'Snowfall': snowfall_list,
            'Sun': sun_list,
            'Wind': wind_list}

    df = pd.DataFrame.from_dict(data, orient='index').transpose()

    # Create HTML table dynamically with cell colors for weather variables
    html_table = """
    <!DOCTYPE html>
    <html>
    <head>
    <style>
      table {
        border-collapse: collapse;
        width: 50%;
        margin: 20px;
      }

      th, td {
        border: 2px solid black;
        padding: 10px;
        text-align: center;
      }

      th {
        background-color: #6b768a;
        color: white;
      }

      tr:nth-child(even) {
        background-color: #f2f2f2;
      }

      tr:hover {
        background-color: #ddd;
      }

      .red {
        background-color: #FF6961;
      }

      .orange {
        background-color: #FAC898;
      }

      .yellow {
        background-color: #fdfd96;
      }

      .lime {
        background-color: #D1FEB8;
      }

      .green {
        background-color: #77DD77;
      }
    </style>
    </head>
    <body>

    <table>
    <caption>Skiing weather variables</caption>
      <tr>{}</tr>
      {}
    </table>

    </body>
    </html>
    """

    # Create table headers dynamically for weather variables
    headers = ('\n'.join(['<th>{}</th>'.format(header)
                          for header in df.columns]))

    # Create table rows dynamically with cell colors according to values
    rows = ''
    for index, row in df.iterrows():
        rows += '<tr>'
        for value in row:
            if value == 0:
                rows += '<td class="red">{}</td>'.format(value)
            elif value == 1:
                rows += '<td class="orange">{}</td>'.format(value)
            elif value == 2:
                rows += '<td class="yellow">{}</td>'.format(value)
            elif value == 3:
                rows += '<td class="lime">{}</td>'.format(value)
            elif value == 4:
                rows += '<td class="green">{}</td>'.format(value)
            else:
                rows += '<td>{}</td>'.format(value)
        rows += '</tr>\n'

    # Using regular expressions, change each of the colored cells with a single
    # integer value to display the correct output data from
    # snowcheck.snow_variables function

    pattern = r'(\b\d\b)'
    rows = re.sub(pattern, get_match_value, rows)

    # Insert headers and rows into the HTML table code
    html_table = html_table.replace("{}", headers, 1).replace("{}", rows, 1)

    return html_table


def geographical_table(lon, lat, ds):
    """
    @authors: Malte Hildebrandt, Joep van Noort

    This function generates an html file containing a table in which certain
    location dependent values provide information about whether or not the
    given location is potentially apt for northern hemisphere winter skiing.

    The function that is used to obtain this information is the mountain_check
    function from the snowcheck package. The table is then dynamically
    generated with value dependent coloring.

    Parameters
    ----------
    lon : user input longitude value in degrees East
    lat : user input latitude value in degrees North
    ds : wrfout dataset

    Returns
    -------
    html_table : html table in string format. A html file is generated in
                 the backgorund.

    """
    # Table 2
    # Create dataset for physical variables

    data2 = {'Minimum Mountain Height': (snowcheck.mountain_check
                                         (lon, lat, ds)[2]),
             'Topographic Height': snowcheck.mountain_check(lon, lat, ds)[3],
             'Snowsure': snowcheck.mountain_check(lon, lat, ds)[4],
             'Skiable Slopes': snowcheck.mountain_check(lon, lat, ds)[5]}

    df2 = pd.DataFrame([data2])

    # Create HTML table for physical attributes
    html_table2 = """
    <!DOCTYPE html>
    <html>
    <head>
    <style>
      table {
        border-collapse: collapse;
        width: 50%;
        margin: 20px;
      }

      th, td {
        border: 2px solid black;
        padding: 10px;
        text-align: center;
      }

      th {
        background-color: #6b768a;
        color: white;
      }

      tr:nth-child(even) {
        background-color: #f2f2f2;
      }

      tr:hover {
        background-color: #ddd;
      }

      .red {
        background-color: #FF6961;
      }

      .green {
        background-color: #77DD77;
      }

    </style>
    </head>
    <body>

    <table>
    <caption>Geographical attributes for skiing in winter</caption>
      <tr>{}</tr>
      {}
    </table>

    </body>
    </html>
    """

    # Create table headers for physical attributes
    headers2 = ('\n'.join(['<th>{}</th>'.format(header)
                           for header in df2.columns]))

    # Create table rows for physical attributes
    rows2 = ''
    for index, row in df2.iterrows():
        rows2 += '<tr>'
        for value in row:
            if value == 'No':
                rows2 += '<td class="red">{}</td>'.format(value)
            elif value == 'Yes':
                rows2 += '<td class="green">{}</td>'.format(value)
            else:
                rows2 += '<td>{}</td>'.format(value)
        rows2 += '</tr>\n'

    # Insert headers and rows into the HTML table code
    html_table2 = (html_table2.replace("{}", headers2, 1).replace("{}",
                                                                  rows2, 1))

    return html_table2


def html_page(html_table, html_table2, filename="tables.html"):
    """
    @ authors: Malte Hildebrandt, Joep van Noort

    This function combines two html tables into one html file that is saved to
    the specified directory.

    Parameters
    ----------
    html_table : input table 1
    html_table2 : input table 2
    filename : name of the HTML file (default is "tables.html")

    Returns
    -------
    None.

    """
    # Create HTML code for the entire page
    html_page = f"""
    <!DOCTYPE html>
    <html>
    <head>
    </head>
    <body>

    <!-- Insert the HTML table code -->
    {html_table}
    {html_table2}

    </body>
    </html>
    """

    output_directory = cfg.output_directory
    # Save the entire HTML page code to a file
    full_path = os.path.join(output_directory, filename)
    with open(full_path, "w") as file:
        file.write(html_page)

    print(f"HTML page with tables saved as {full_path}")


# ################### Test run ####################

# test values
ds = xr.open_dataset(cfg.wrfout)
lon = 11
lat = 47
time = 24
counter = 1
timestep = 0

html_table = weather_table(lon, lat, ds)
html_table2 = geographical_table(lon, lat, ds)
html_page(html_table, html_table2)
