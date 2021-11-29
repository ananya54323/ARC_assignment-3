# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 21:25:40 2021

Name: Ananya Balaji
ID: 17280049

"""

#!/usr/bin/python

import os, sys
import json
import numpy as np
import re
from pathlib import Path
from matplotlib import colors
import matplotlib.pyplot as plt


def solve_0a938d79(x): #this is a flash fill example where one element is colored and that decides the colour of the elements in that column
    def row_flash_fill(grid): #filling rows and checking the similarity
        n_colors = [] #getting the distinct colors present in the task
        color_index = [] #non black colours
        for i in range(len(grid)):
            n_colors.append(len(np.unique(grid[i])))
        for i in range(len(n_colors)):
            if n_colors[i] == max(n_colors):
                color_index.append(i)
        diff = color_index[1]-color_index[0] #finding the delta to assess the difference in the colours
        #flash filling the first row
        for i in range(color_index[0], len(n_colors), 2*diff):
            for j in range(len(grid[0])):
                grid[i][j] = np.unique(grid[color_index[0]])[-1]
        #flash fill second row    
        for i in range(color_index[1], len(n_colors), 2*diff):
            for j in range(len(grid[0])):
                    grid[i][j] = np.unique(grid[color_index[1]])[-1]
        return grid
    def col_flash_fill(grid):
        n_colors = [] #number of distinct colours in columns
        color_index=[] #non black columns
        column = np.array(grid).T
        for j in range(len(column)):
            n_colors.append(len(np.unique(column[j])))
        for i in range(len(n_colors)):
            if n_colors[i] == max(n_colors):
                color_index.append(i)     
        diff = color_index[1] - color_index[0]
        #first column flash fill
        for j in range(color_index[0],len(n_colors),2*diff):
            for i in range(len(grid)):
                grid[i][j] = np.unique(column[color_index[0]])[-1]
            #second column flash fill
            for j in range(color_index[1], len(n_colors),2*diff):
                for i in range(len(grid)):
                    grid[i][j] = np.unique(column[color_index[1]])[-1]
        return grid
    #print output
    if len(x)>len(x[0]):
        row_flash_fill(x)
    elif len(x)<len(x[0]):
        col_flash_fill(x)
    return x

def solve_0b148d64(x):#this is a cropping task where one entity is cropped from the 4 diagrams in the tasks.
    def min_crop(grid):
        x = np.bincount(grid.flatten(),minlength=10) #counting number of occurences
        y = int(np.where(x==np.min(x[np.nonzero(x)]))[0])
        coord = np.argwhere(grid==y)#indicesof non zero elements, grouped by element
        xmin,ymin= coord.min(axis=0)
        xmax,ymax = coord.max(axis =0)
        return grid[xmin:xmax+1, ymin:ymax+1]
    grid = min_crop(x)
    return grid 

def main():
    # Find all the functions defined in this file whose names are
    # like solve_abcd1234(), and run them.

    # regex to match solve_* functions and extract task IDs
    p = r"solve_([a-f0-9]{8})" 
    tasks_solvers = []
    # globals() gives a dict containing all global names (variables
    # and functions), as name: value pairs.
    for name in globals(): 
        m = re.match(p, name)
        if m:
            # if the name fits the pattern eg solve_abcd1234
            ID = m.group(1) # just the task ID
            solve_fn = globals()[name] # the fn itself
            tasks_solvers.append((ID, solve_fn))

    for ID, solve_fn in tasks_solvers:
        # for each task, read the data and call test()
        directory = os.path.join("C:\\","users","ananya","arc","data", "training")
        json_filename = os.path.join(directory, ID + ".json")
        data = read_ARC_JSON(json_filename)
        test(ID, solve_fn, data)
    
def read_ARC_JSON(filepath):
    """Given a filepath, read in the ARC task data which is in JSON
    format. Extract the train/test input/output pairs of
    grids. Convert each grid to np.array and return train_input,
    train_output, test_input, test_output."""
    
    # Open the JSON file and load it 
    data = json.load(open(filepath))

    # Extract the train/test input/output grids. Each grid will be a
    # list of lists of ints. We convert to Numpy.
    train_input = [np.array(data['train'][i]['input']) for i in range(len(data['train']))]
    train_output = [np.array(data['train'][i]['output']) for i in range(len(data['train']))]
    test_input = [np.array(data['test'][i]['input']) for i in range(len(data['test']))]
    test_output = [np.array(data['test'][i]['output']) for i in range(len(data['test']))]

    return (train_input, train_output, test_input, test_output)


def test(taskID, solve, data):
    """Given a task ID, call the given solve() function on every
    example in the task data."""
    print(taskID)
    train_input, train_output, test_input, test_output = data
    print("Training grids")
    for x, y in zip(train_input, train_output):
        yhat = solve(x)
        show_result(x, y, yhat)
    print("Test grids")
    for x, y in zip(test_input, test_output):
        yhat = solve(x)
        show_result(x, y, yhat)

        
def show_result(x, y, yhat):
    print("Input")
    print(x)
    print("Correct output")
    print(y)
    print("Our output")
    print(yhat)
    print("Correct?")
    if y.shape != yhat.shape:
        print(f"False. Incorrect shape: {y.shape} v {yhat.shape}")
    else:
        print(np.all(y == yhat))


if __name__ == "__main__": main()
