from raresim.common.exceptions import RaresimException, IllegalArgumentException
from numpy import ndarray
import random
import os
import gzip
import timeit
import numpy as np
import numba as nb
import array
import sys

class SparseMatrix:
    """
    A sparse matrix is a matrix of 0s and 1s in which most of the values are 0s.
    It stores the data efficiently by simply storing the indices of each 1 rather than storing all the 1s and 0s
    of the matrix. The __data field is a list of "rows" where each "row" is a list the indices where a 1 is present
    in that row.
    """

    def __init__(self, cols=0):
        self.__cols = cols
        self.__rows = 0
        self.__data = []

    def set_col_count(self, cols: int) -> None:
        """
        Only to be used by binary reader to set the number of columns when creating SparseMatrix object
        @param cols: number of cols
        @return: None
        """
        self.__cols = cols

    def get(self, row: int, col: int) -> int:
        """
        Gets the value in the matrix at a requested position
        @param row: index of row
        @param col: index of column
        @return: value at requested position
        """
        if row > self.__rows or col > self.__cols:
            raise RaresimException(
                f"Attempted to get value at {row},{col} when the bounds of the matrix were {self.__rows},{self.__cols}")

        arr = self.__data[row]
        # Use binary search to search for the presence of the index in O(log(n)) time
        low = 0
        high = len(arr) - 1
        while low <= high:
            mid = (high + low) // 2
            if arr[mid] < col:
                low = mid + 1
            elif arr[mid] > col:
                high = mid - 1
            else:
                # Return 1 if index is present
                return 1
        # If we reach here, then the index was not present so the value at the requested position is a 0
        return 0

    def get_row(self, row: int) -> list:
        """
        Gets the full row at the requested index
        @param row: index of row to get
        @return: the full row at the requested index
        """
        if row > self.__rows:
            raise RaresimException(f"Attempted to get row {row}, but there are only {self.__rows} rows")
        ret = [0] * self.__cols
        for i in self.__data[row]:
            ret[i] = 1
        return ret

    def get_row_raw(self, row: int) -> list:
        """
        Gets the indices of the ones in a row
        @param row: index of row to get
        @return: list of the indices of the ones in a row
        """
        return self.__data[row]

    def add(self, row: int, col: int) -> None:
        """
        Adds a 1 to the matrix at the given location.
        @param row: row to insert the value
        @param col: column to insert the one at
        @return:
        """
        if col > self.__cols:
            raise Exception(f"Attempted to insert at column {col}, but the matrix only has {self.__cols} columns")
        while row > len(self.__data):
            self.__data.append(())
            self.__rows += 1
        temp = self.__data[row]
        temp.append(col)
        temp.sort()
        temp = list(set(temp))
        self.__data[row] = temp

    def remove(self, row: int, col: int) -> None:
        """
        Places a zero at the given position
        @param row: row index
        @param col: column index
        @return: None
        """
        self.__data[row].remove(col)


    def add_row(self, val: list) -> None:
        """
        Sets row at given index to the provided list. Only to be used by SparseMatrixReader
        @param val: list of indices representing locations of 1s in the new row
        @return: None
        """
        self.__data.append(val)
        self.__rows += 1

    def remove_row(self, row: int) -> None:
        """
        Removes the row at the given index
        @param row: index of the row to remove
        @return: None
        """
        self.__data.pop(row)
        self.__rows -= 1

    def num_rows(self) -> int:
        """
        @return: the number of rows in the matrix
        """
        return self.__rows

    def num_cols(self) -> int:
        """
        @return: the number of columns in the matrix
        """
        return self.__cols

    def row_num(self, row: int) -> int:
        """
        @param row: index of row
        @return: the number of 1s in the requested row
        """
        if row > self.__rows - 1:
            raise RaresimException(f"Attempted to access row {row} but there were only{self.__rows} in the matrix")
        return len(self.__data[row])

    def prune_row(self, row: int, num_prune: int) -> None:
        """
        Randomly prunes n 1s from the given row
        @param row: index of the row to prune
        @param num_prune: how many 1s to prune from the row
        @return: None
        """
        if row > self.__rows:
            return

        num_keep = len(self.__data[row]) - num_prune
        keep_ids = self.__reservoir_sample(num_keep, len(self.__data[row]))
        keep_ids.sort()
        ret = []
        for i in keep_ids:
            ret.append(self.__data[row][i])

        self.__data[row] = ret

    @staticmethod
    def __reservoir_sample(k: int, n: int) -> list:
        """
        @param k: desired size of sample
        @param n: max size of any element + 1
        @return: random sample of k elements in the range [0...n-1]
        """
        stream = [i for i in range(n)]
        reservoir = [0] * k
        i = 0
        for i in range(k):
            reservoir[i] = stream[i]

        while i < n:
            j = random.randrange(i + 1)
            if j < k:
                reservoir[j] = stream[i]
            i += 1
        return reservoir

class SparseMatrixReader:
    def __init__(self):
        pass

    def loadSparseMatrix(self, filepath: str) -> SparseMatrix:
        """
        Loads a sparse matrix from the given file
        @param filepath: path to file to read
        @return: loaded sparse matrix
        @throws IllegalArgumentException: if the file does not exist
        """
        if not os.path.isfile(filepath):
            raise IllegalArgumentException(f"No such file exists: {filepath}")

        if filepath[-3:] == '.sm':
            ret = self.__loadCompressed(filepath)
        elif filepath[-3:] == '.gz':
            ret = self.__loadZipped(filepath)
        else:
            ret = self.__loadUncompressed(filepath)
        return ret

    def __loadZipped(self, filepath: str) -> SparseMatrix:
        """
        Loads a sparse matrix from the given g-zipped file. The resulting sparse matrix will have the same number of rows and columns
        as the input file.
        @param filepath: path to file to read
        @return: loaded sparse matrix
        """
        matrix = None
        for line in gzip.open(filepath, "rt"):
            if line is None or line.strip() == '':
                break
            nums = self.compute(np.frombuffer(line[::2].encode('ascii'), np.uint8))
            if matrix is None:
                matrix = SparseMatrix(len(nums))
            row_to_add = self.__getSparseRow(nums)
            matrix.add_row(row_to_add.tolist())

        return matrix

    def __loadCompressed(self, filepath: str) -> SparseMatrix:
        """
        Loads a sparse matrix from the given binary encoded file. See readme for details on encoding.
        Values of -1 (0xFFFFFFFF) are considered row delimiters.
        @param filepath: path to file to read
        @return: loaded sparse matrix
        """
        with open(filepath, "rb") as f:
            data = f.read(4)
            matrix = SparseMatrix(int.from_bytes(data, "little"))
            row = []
            data = f.read(4)
            while data:
                if self.__toSigned32(int.from_bytes(data, "little")) == -1:
                    matrix.add_row(np.fromiter(row, dtype=int).tolist())
                    row = []
                else:
                    row.append(int.from_bytes(data, "little"))
                data = f.read(4)
        return matrix

    def __loadUncompressed(self, filepath: str) -> SparseMatrix:
        """
        Loads a sparse matrix from the given uncompressed human-readable file.
        The resulting sparse matrix will have the same number of rows and columns
        as the input file.
        @param filepath: path to file to read
        @return: loaded sparse matrix
        """
        matrix = None
        for line in open(filepath, "r"):
            if line is None or line.strip() == '':
                break
            nums = self.compute(np.frombuffer(line[::2].encode('ascii'), np.uint8))
            if matrix is None:
                matrix = SparseMatrix(len(nums))
            row_to_add = self.__getSparseRow(nums)
            matrix.add_row(row_to_add.tolist())

        return matrix

    @staticmethod
    def __getSparseRow(nums: np.ndarray) -> ndarray:
        return np.where(nums == 1)[0]

    @staticmethod
    def __toSigned32(n):
        n = n & 0xffffffff
        return n | (-(n & 0x80000000))

    @staticmethod
    @nb.njit(nb.int32[::1](nb.types.Array(nb.uint8, 1, 'C', readonly=True)))
    def compute(arr):
        """
        This method is not of my own design, but is the fastest possible way that I know of (without writing my own
        C-based extension) to convert a long delimited string into a list of ints. I found this solution at
        https://stackoverflow.com/questions/74873414/the-fastest-way-possible-to-split-a-long-string

        Due to this method being 'pre-compiled' by the numba just-in-time compiler, it cannot be debugged unless the
        @nb.njit decorator line is commented out. This will cause a very noticeable performance degradation, but will
        allow for the placement of breakpoints in this method for debugging.
        """
        count = len(arr)
        res = np.empty(count, np.int32)
        base = ord('0')
        val = 0
        cur = 0
        for c in arr:
            val = (val * 10) + c - base
            res[cur] = val
            cur += 1
            val = 0
        return res

class SparseMatrixWriter:
    def __init__(self):
        pass

    def writeToHapsFile(self, sparseMatrix: SparseMatrix, filename: str, compression="gz") -> None:
        """
        Writes the given sparse matrix to a file with the given name and compression method.

        @param sparseMatrix: input matrix
        @param filename: output file
        @param compression: compression method. Default is "gz". Can be "gz" for g-zipped, "sm" for binary encoded, or "" for uncompressed.
        @return: None
        """
        if compression == "gz":
            self.__writeZipped(sparseMatrix, filename)
        elif compression == "sm":
            self.__writeCompressed(sparseMatrix, filename)
        else:
            self.__writeUncompressed(sparseMatrix, filename)
        sys.stdout.write("\r[%-20s] %d%%" % ('='* 20, 100))
        print()

    @staticmethod
    def __writeZipped(sparseMatrix: SparseMatrix, filename: str):
        """
        Writes the sparse matrix to a g-zipped format. Unzipping will yield a human-readable file
        @param sparseMatrix: input matrix
        @param filename: output file
        @return: None
        """
        with gzip.open(filename, "wb") as f:
            for i in range(sparseMatrix.num_rows()):
                row = ["0"]*sparseMatrix.num_cols()
                for j in sparseMatrix.get_row_raw(i):
                    row[j] = "1"
                line = " ".join(row) + "\n"
                f.write(line.encode())

                sys.stdout.write("\r[%-20s] %d%%" % ('='* int((i / sparseMatrix.num_rows()) * 20), i/sparseMatrix.num_rows()*100))

    @staticmethod
    def __writeUncompressed(sparseMatrix: SparseMatrix, filename: str):
        """
        Writes the sparse matrix to an uncompressed human-readable format
        @param sparseMatrix: input matrix
        @param filename: output file
        @return: None
        """
        step = int(sparseMatrix.num_rows() / 10)
        with open(filename, "w") as f:
            for i in range(sparseMatrix.num_rows()):
                row = ["0"] * sparseMatrix.num_cols()
                for j in sparseMatrix.get_row_raw(i):
                    row[j] = "1"
                line = " ".join(row) + "\n"
                f.write(line)
                if i % step == 0:
                    print('.', end='', flush=True)


    @staticmethod
    def __writeCompressed(sparseMatrix: SparseMatrix, filename: str):
        """
        Writes the sparse matrix to a binary encoded file. Values of -1 (0xFFFFFFFF) are considered row delimiters.
        @param sparseMatrix: input matrix
        @param filename: output file
        @return: None
        """
        step = int(sparseMatrix.num_rows() / 10)
        with open(filename, "wb") as f:
            f.write(int.to_bytes(sparseMatrix.num_cols(), 4, "little"))
            for i in range(sparseMatrix.num_rows()):
                row = sparseMatrix.get_row_raw(i)
                data = array.array("i", row + [-1])
                f.write(data.tobytes())
                if i % step == 0:
                    print('.', end='', flush=True)
