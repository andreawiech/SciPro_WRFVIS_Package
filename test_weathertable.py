import xarray as xr
import unittest
import re

from wrfvis.tables import weather_table
from wrfvis import cfg


class TestWeatherTableFunction(unittest.TestCase):
    """
    @author: Joep van Noort

    This class contains 3 tests to check that the weather_table function
    creates a correctly looking html. The color-coded values are tested for
    correctness for the snowfall, sunshine and wind columns. Test location
    is in Tyrol, Austria.
    """
    # test values
    ds = xr.open_dataset(cfg.wrfout)
    lon = 11
    lat = 47
    # Call the function
    html_table = weather_table(lon, lat, ds)

    def test_weather_table_snow(self):
        """
        @author: Joep van Noort

        This is a test to check if snowfall values displayed in the
        weather_tables html table display in the correct colored cell.

        """
        # Use Regex to find all cells containing snowfall
        snowfall_sequences = re.findall(
            r'(<td class="\w+">\d+ cm</td>)', self.html_table)

        for sequence in snowfall_sequences:
            # Extract the numerical value from the string
            value = float(sequence.split('>')[1].split()[0])

            # Get the corresponding class based on the value
            if 0 < value < 5:
                expected_class = 'orange'
            elif 5 <= value < 15:
                expected_class = 'yellow'
            elif 15 <= value < 30:
                expected_class = 'lime'
            elif 30 <= value:
                expected_class = 'green'
            else:
                expected_class = 'red'

            # Extract the class from the sequence
            actual_class = sequence.split('"')[1]

            # Remove the test-specific part of the class
            actual_class = actual_class.split()[0]

            # Compare the expected and actual classes
            assert expected_class == actual_class, (
                f"Failed for value {value}: "
                f"Expected {expected_class}, but got {actual_class}"
                )

    def test_weather_table_sunshine(self):
        """
        @author: Joep van Noort

        This is a test to check if sunshine values displayed in the
        weather_tables html table display in the correct colored cell.

        """

        # Define the accepted sequences
        accepted_sequences = [
            '<td class="green">Sunny</td>',
            '<td class="lime">Mostly Sunny</td>',
            '<td class="yellow">Somewhat Sunny</td>',
            '<td class="orange">Mostly Cloudy</td>',
            '<td class="red">Overcast</td>'
        ]

        # Check if at least one accepted sequence is present
        self.assertTrue(any(sequence in self.html_table
                            for sequence in accepted_sequences),
                        "At least one accepted sequence should "
                        "be present in the HTML table.")

        # Check that no unexpected sequences are present
        unexpected_sequences = [
            f'<td class="{weather_class}">{condition}</td>'
            for weather_class in ['green', 'lime', 'yellow', 'orange', 'red']
            for condition in ['Sunny', 'Mostly Sunny', 'Somewhat Sunny',
                              'Mostly Cloudy', 'Overcast']
            if f'<td class="{weather_class}">{condition}</td>'
            not in accepted_sequences
        ]
        for sequence in unexpected_sequences:
            self.assertNotIn(
                sequence, self.html_table,
                f"Unexpected sequence found: {sequence}")

    def test_weather_table_wind(self):
        """
        @author: Joep van Noort

        This is a test to check if windspeed values displayed in the
        weather_tables html table display in the correct colored cell.

        """
        # Use Regex to find all cells containing windspeed
        windspeed_sequences = re.findall(
            r'(<td class="\w+">\d+.\d+ m/s</td>)', self.html_table)

        for sequence in windspeed_sequences:
            # Extract the numerical value from the string
            value = float(sequence.split('>')[1].split()[0])

            # Get the corresponding class based on the value
            if 0.00 <= value <= 1.00:
                expected_class = 'green'
            elif 1.00 < value < 3.00:
                expected_class = 'lime'
            elif 3.00 <= value < 8.00:
                expected_class = 'yellow'
            elif 8.00 <= value < 15.00:
                expected_class = 'orange'
            else:
                expected_class = 'red'

            # Extract the class from the sequence
            actual_class = sequence.split('"')[1]

            # Remove the test-specific part of the class
            actual_class = actual_class.split()[0]

            # Compare the expected and actual classes
            assert expected_class == actual_class, (
                f"Failed for value {value}: "
                f"Expected {expected_class}, but got {actual_class}"
                )


if __name__ == '__main__':
    unittest.main()
