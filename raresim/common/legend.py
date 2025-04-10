import os
from raresim.common.exceptions import IllegalArgumentException
import warnings

class Legend:
    supported_columns = ["id", "position", "a0", "a1", "prob", "protected", "fun", "AC", "exonic", "gene"]

    def __init__(self, header: list):
        self.__header = header
        self.__rows = []

    def get_header(self) -> list:
        """
        Returns the column names of the legend

        The column names are the same as the header of the legend file that was loaded to create this object.

        :return: a list of strings that are the column names of the legend
        """
        return self.__header

    def row_count(self) -> int:
        """
        Returns the number of rows in the legend.

        :return: the number of rows in the legend
        """
        return len(self.__rows)

    def add_row(self, row: list) -> None:
        """
        Adds a row to the legend.

        The row must have the same length as the header of the legend, and each element of the row must match the type of the corresponding element in the header.

        :param row: the row to add
        :return: None
        """
        self.__rows.append(row)

    def remove_row(self, index: int) -> None:
        """
        Removes a row from the legend.

        The row at the given index will be removed from the legend.

        :param index: the index of the row to remove
        :return: None
        """
        self.__rows.pop(index)

    def get_row(self, index: int) -> dict:
        """
        Returns the row at the given index as a dictionary.

        The keys of the dictionary are the column names of the legend, and the values are the values of the corresponding elements in the row.

        :param index: the index of the row to get
        :return: the row at the given index as a dictionary
        """
        return dict(zip(self.__header, self.__rows[index]))

    def __getitem__(self, index: int):
        """
        Returns the row at the given index as a dictionary.

        The keys of the dictionary are the column names of the legend, and the values are the values of the corresponding elements in the row.

        :param index: the index of the row to get
        :return: the row at the given index as a dictionary
        """
        return self.get_row(index)

class LegendReaderWriter:
    def __init__(self):
        pass

    @staticmethod
    def load_legend(filepath: str) -> Legend:
        """
        Loads a legend file from the given filepath and returns a Legend object.

        The legend file should have a header line with column names, and each subsequent line
        should contain a row of data. The column names should be one of the following:
        "id", "position", "a0", "a1", "prob", "protected", or "fun".

        If the file does not exist, an IllegalArgumentException is raised.
        If a column name is not recognized, a RaresimException is raised.

        :param filepath: path to the legend file
        :return: a Legend object
        """
        if not os.path.isfile(filepath):
            raise IllegalArgumentException(f"No such file exists: {filepath}")

        with open(filepath, "r") as f:
            line = f.readline()
            header = line.rstrip().split()
            for key in header:
                if key not in Legend.supported_columns:
                    warnings.warn(f"Legend column '{key}' is not supported. Supported keys are {Legend.supported_columns}")
            legend = Legend(header)

            line = f.readline()
            while line and line.strip() != "\n" and line.strip() != '':
                row = line.rstrip().split('\t')
                legend.add_row(row)
                line = f.readline()
        return legend

    @staticmethod
    def write_legend(legend: Legend, filepath: str) -> None:
        """
        Writes the given Legend object to a file at the given filepath.

        The file will have a header line with column names, and each subsequent line
        will contain a row of data. The column names will be the same as the Legend
        object's header.

        If the file cannot be opened for writing, an IOError is raised.

        :param legend: the Legend object to write
        :param filepath: the path to the file to write
        """
        with open(filepath, "w") as f:
            header_string = "\t".join(legend.get_header()) + "\n"
            f.write(header_string)
            for i in range(legend.row_count()):
                line = "\t".join(legend.get_row(i)) + "\n"
                f.write(line)
