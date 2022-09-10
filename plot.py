import matplotlib.pyplot as plt

# line 1 points
x1 = [10,20,30,40,50]
y1 = [324.9,1585.34,2106.76,2642.96,3108.18]
# plotting the line 1 points
plt.plot(x1, y1, label="WL", marker='s', linestyle='dashed',markersize=4)

# line 2 points
x2 = [10,20,30,40,50]
y2 = [231.54,1468.54,2124.14,2485.74,3254.46]
# plotting the line 2 points
plt.plot(x2, y2, label="WS", marker='o', markersize=4)

# naming the x axis
plt.xlabel('Problem size, n')
# naming the y axis
plt.ylabel('Average Run Time (s)')
# giving a title to my graph
plt.title('Average run time on dataset 2 instances')

# show a legend on the plot
plt.legend()

# function to show the plot
plt.show()