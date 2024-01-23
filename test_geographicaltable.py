import xarray as xr
import unittest
import re

from wrfvis.tables import geographical_table
from wrfvis import cfg


class TestGeographicalTableFunction(unittest.TestCase):
    """
    @author: Joep van Noort

    This class contains 2 tests to check that the geographical_table function
    creates a correctly looking html. The color-coded values are tested for
    correctness for the Snowsure and Skiable Slope cells, topographic height
    is checked to be under a maximum and minimum mountain height is checked
    to be positive. Test location is in Tyrol, Austria.
    """
    # test values
    ds = xr.open_dataset(cfg.wrfout)
    lon = 11
    lat = 47
    # Call the function
    html_table = geographical_table(lon, lat, ds)

    def test_geographical_table_snowslope(self):
        """
        @author: Joep van Noort

        This is a test to check if snowsure and slope values displayed in the
        geographical_tables html table display in the correct colored cell.

        """
        # Use Regex to find all cells containing snowfall
        snowslope_sequences = re.findall(
            r'(<td class="\w+">\w+</td>)', self.html_table)

        for sequence in snowslope_sequences:
            # Extract the Yes-No value from the string
            value = sequence.split('>')[1].split('<')[0]

            # Get the corresponding class based on the value
            if value == "Yes":
                expected_class = 'green'
            elif value == "No":
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

    def test_geographical_table_height(self):
        """
        @author: Joep van Noort

        This is a test to check if mountain height values displayed in the
        geographical_tables html table are realistic.

        """

        # Use Regex to find all cells containing snowfall
        height_sequences = re.findall(
            r'(<td>\d,\d+ meters</td>)', self.html_table)

        for sequence in height_sequences:
            # Extract the height value from the string as an integer
            value = int(''.join(sequence.split('>')[1].split()[0].split(',')))

            # Check if value is positive and under 4810 (Height of Mont Blanc)
            assert value >= 0, (f"Failed for sequence {sequence}: "
                                f"Expected positive height, but got {value}")
            assert value <= 4810, (f"Failed for sequence {sequence}: "
                                   f"There are no such heights in W-Europe")


if __name__ == '__main__':
    unittest.main()
