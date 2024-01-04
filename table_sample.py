import pandas as pd
import xarray as xr
import re

from wrfvis import cfg, snowcheck


def match_counter(match):
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
    Time_list = []
    Snowcover_list = []
    Snowfall_list = []
    Sun_list = []
    Wind_list = []

    for i in timestep:
        date = ds['XTIME'][i].values
        time_value = date.astype(str)[5:10] + ' ' + date.astype(str)[11:16]
        Time_list.append(time_value)
        Snowcover_list.append(snowcheck.snow_variables(lon, lat, ds, i)[0])
        Snowfall_list.append(snowcheck.snow_variables(lon, lat, ds, i)[2])
        Sun_list.append(snowcheck.snow_variables(lon, lat, ds, i)[4])
        Wind_list.append(snowcheck.snow_variables(lon, lat, ds, i)[6])

    data = {'Time': Time_list,
            'Snowcover': Snowcover_list,
            'Snowfall': Snowfall_list,
            'Sun': Sun_list,
            'Wind': Wind_list}

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
    rows = re.sub(pattern, match_counter, rows)

    # Insert headers and rows into the HTML table code
    html_table = html_table.replace("{}", headers, 1).replace("{}", rows, 1)

    # Save HTML table code to a file
    table_filename = "colored_cells_table.html"
    with open(table_filename, "w") as file:
        file.write(html_table)

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

    # Save HTML table code to a file
    table_filename2 = "physical_attributes_table.html"
    with open(table_filename2, "w") as file:
        file.write(html_table2)

    return html_table2


def html_page(html_table, html_table2):
    """
    @ authors: Malte Hildebrandt, Joep van Noort

    This function combines two html tables into one html file that is saved to
    the local directory.

    Parameters
    ----------
    html_table : input table 1
    html_table2 : input table 2

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

    # Save the entire HTML page code to a file
    html_filename = "combined_page.html"
    with open(html_filename, "w") as file:
        file.write(html_page)

    print(f"HTML page with tables saved as {html_filename}")


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
