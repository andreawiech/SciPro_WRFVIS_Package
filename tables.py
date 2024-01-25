"""
This script provides functions to generate HTML tables containing weather and
geographical information for skiing conditions based on WRF (Weather Research
and Forecasting) model output.

Authors: Malte Hildebrandt, Joep van Noort

Dependencies:
- pandas
- wrfvis package (imported as snowcheck)
- Beautiful Soup 4

Note: Ensure that the necessary dependencies are installed before
running the script.
"""
import pandas as pd
from bs4 import BeautifulSoup
import os

from wrfvis import snowcheck, cfg


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

    # create timeseries list of weather variable values
    timestep = range(36)
    replacement_list = []

    for i in timestep:
        replacement_list.append(snowcheck.snow_variables(lon, lat, ds, i)[1])
        replacement_list.append(snowcheck.snow_variables(lon, lat, ds, i)[3])
        replacement_list.append(snowcheck.snow_variables(lon, lat, ds, i)[5])
        replacement_list.append(snowcheck.snow_variables(lon, lat, ds, i)[7])

    soup = BeautifulSoup(rows, 'html.parser')

    # Find all <td> tags with class="red", "green", or "yellow"
    td_tags = soup.find_all('td', class_=['red', 'orange',
                                          'lime', 'green', 'yellow'])

    # Iterate through the td_tags and replace the text content with elements
    # from the replacement_list
    for i, td_tag in enumerate(td_tags):
        td_tag.string.replace_with(str(replacement_list[i % len(replacement_list)]))

    # Assign the modified HTML string
    rows = str(soup)

    # Load HTML template from file
    with open(cfg.html_weather_template, 'r') as template_file:
        html_template = template_file.read()

    # Insert headers and rows into the HTML template code
    html_table = html_template.replace("{}", headers, 1).replace("{}", rows, 1)

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

    # Load HTML template from file
    with open(cfg.html_geographic_template, 'r') as template_file:
        html_template2 = template_file.read()

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
    html_table2 = (html_template2.replace("{}", headers2, 1)
                    .replace("{}", rows2, 1))

    return html_table2


def html_page(html_table, html_table2, filename="snowcheck.html", directory=None):
    """
    @ authors: Malte Hildebrandt, Joep van Noort

    This function combines two html tables into one html file that is saved to
    the specified or default local directory.

    Parameters
    ----------
    html_table : input table 1
    html_table2 : input table 2
    filename : str, optional
        Name of the HTML file. Default is "snowcheck.html".
    directory : str, optional
        Directory where the HTML file will be saved. If not provided, a default directory may be used.

    Returns
    -------
    None.
    """

    # Create directory for the plot if not provided
    if directory is None:
        directory = cfg.output_directory

    # Load HTML template from file
    with open(cfg.html_twotable_template, 'r') as template_file:
        html_template_page = template_file.read()

    # Save the entire HTML page code to a file in the specified directory
    full_path = os.path.join(directory, filename)
    with open(full_path, "w") as file:
        file.write(html_template_page.format(html_table, html_table2))
