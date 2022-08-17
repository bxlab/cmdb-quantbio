---
title: Numpy Cheat Sheet
layout: default
---

<!-- Make columns even -->
<style>
  td:first-child { width: 20% }
</style>

## Numpy Cheat Sheet

#### ARRAYS
##### An ordered collection of elements that can be modified (mutable) and retrieved (indexed). This collection may be multi-dimensional.

#### Array creation

|:-----------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `myArray = numpy.array(['abc', 'bca', 'cab'])` | Create a new one-dimensional array called 'myArray' with string elements. |
| `my2DArray = numpy.array([[0, 15, 4],[8, 12, 0]])` | Create a new two-dimensional array called `my2DArray` with integer elements |
| `newArr = numpy.array(variable)`             | Convert variable into a numpy array |
| `newEmptyArr = numpy.empty((shape1, shape2), dtype=np.float32)` | Create a new two-dimensional array that is empty but will store floats shape1*shape2 floats where shape1 is number of rows |
| `newZerosArr = numpy.zeros((shape1, shape2))` | Create a new two-dimensional array that contnains only zeros, where shape1 is the number of rows |
| `newFullArr = numpy.full((shape1, shape2), value)` | Create a new two-dimensional array that contains some value `value`, like `numpy.nan` or `42`; shape1 is the number of rows | 
{:.table.table-striped}                                                                                                                                                         

#### Subsetting and data-selection

:-----------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `retrieved = myArray[index]`                       | Retrieve an item from the one-dimensional `myArray` at the given index. You can change the item at the given index by setting it equal to a new value. An index of -1 will return the last item in the list. **Indexing begins at 0** |
| `retrieved = my2DArray[index]`                     | Retrieve a row of elements at the given index (dtype is another numpy array) from the two-dimensional `my2DArray`. **Indexing begins at 0** |
| `retreived = my2DArray[:, index]`                  | Retrieve a column of elements at the given index (dtype is another numpy array) from the two-dimensionnal `my2DArray`. **Indexing begins at 0** |
| `retrieved = my2DArray[index1, index2]`  | Retrieve a single element at the given indices from the two-dimensional `my2DArray`. **Indexing begins at 0** |
| `slicedArr = myArray[start:stop]`     | Create a new array from a slice of myArray. Either start or stop can be ommitted (e.g. myArray[:3] will give you a slice of myArray up to index 3). Indexing again begins at 0                                       |
{:.table.table-striped}

#### Array attributes

|:-----------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `myArray.ndim`             | Return the number of dimensions of `myArray` |
| `myArray.shape`            | Return the shape (or number of elements) of each dimension of `myArray` as a tuple. (e.g. `myArray.shape[0]` will give you the shape of the first dimension) | 
{:.table.table-striped}

#### Operations and inspecting data

|:-----------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `my2DArray += value`                  | Add an int or float value to every element in the array `my2DArray` |
| `min_value = numpy.amin(my2DArray, axis=0)` | Return the minimum value for each column in the 2D array |
| `min_value2 = numpy.amin(my2DArray, axis=1)` | Return the minimum value for each row in the 2D array |
| `min_value3 = numpy.amin(my2DArray)`  | Return the minimum value in the whole 2D array |
| `max_value = numpy.amax(my2DArray, axis=0)` | Return the maximum value for each column in the 2D array |
| `max_value2 = numpy.amax(my2DArray, axis=1)` | Return the maximum value for each row in the 2D array |
| `max_value3 = numpy.amax(my2DArray)`  | Return the maximum value in the whole 2D array |
| `min_val = numpy.minimum(my2DArray, value)` | Returns the minimum value element-wise, comparing `my2DArray` and the given value, or another array of the same shape. *If any NaNs are present, that is the default minimum* |
| `max_val = numpy.maximum(my2DArray, value)` | Returns the maximum value element-wise, comparing `my2DArray` and the given value, or another array of the same shape. *If any NaNs are present, that is the default maximum* |
| `indices_true = numpy.where(my2DArray > value)` or `indices_true = numpy.where(my2DArray == value2)` | Returns a tuple with arrays of the indices where the condition (`> value`, `== value2`, etc.) is true. If you want the array of row indices, you would need to do `indices_true[0]`. Likewise, for the columns, `indices_true[1]`. Will return empty arrays within the tuple if no elements meet the condition. Will also return a tuple, even for 1D arrays. |              
{:.table.table-striped}

