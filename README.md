# About 
This repository presents two hybrid Genetic Algorithms to solve a combination of production scheduling and bin packing problems in the printing industry. 
To learn about the problem and the algorithms, refer to the following paper

[M.Mostajabdaveh, S.Salman, N.Tahmasbi, Two dimensional guillotine cutting stock and scheduling problem in printing industry, Computers & Operations Research, 2022](https://doi.org/10.1016/j.cor.2022.106014)

#### Please cite the above paper if you are using the algorithm or data.

# Input data
The raw dataset from an actual printing company is available at [Mendely Data](10.17632/bxh46tps75.5) . 

From this dataset, we generated 335 problem instances which can be found in `Data/Input` directory. 

## Input file format
The input files are binary files generated by pickle python package storing an object of Input class (please see the Input class in `Input.py` file)

## Generate
You can generate more problem instances using `Raw Data/Big_Data.xlsx` in  [Mendely Data](10.17632/bxh46tps75.5) and the `Data/Data_Gen.py` . 

# Usage
You can use either RGA or UGA to solve the problem. Generally, RGA is faster, and UGA provides a better solution. Refer to (the paper)[] for a detailed explanation and comparison. 

To use the algorithm, you can simply run the following command

```
python3 runner.py <RGA,UGA> <inputFile> --rep <int> --draw
```

* First input is the type of algorithm which should be either `RGA` or `UGA`
* Second input is the input file name located at Data/Input
* The optional input `--rep <int>` is the number of times we run the algorithm on the input data (default value is 1). In the case of rep>=2, the output result is the average of the best solution costs and runtime. 
* The optional flag `--draw`, if given, the algorithm will print detailed information about the best solution and draw the layout of the bin. 


## Parameters
The algorithm parameters are tuned using a subset of generated instances. However, if you wish to change their values, refer to the `Parameters` class (line 15 runner.py file). 

<!-- CONTRIBUTING -->
## Contributing

Contributions make the open source community an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion to improve this, please fork the repo and create a pull request. You can also open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request





<!-- LICENSE -->
## License

It is distributed under the  GPL-3.0 License. See `LICENSE.txt` for more information.





<!-- CONTACT -->
## Contact

Your Name - [mahdi.ms86@gmail.com](mahdi.ms86@gmail.com)

Project Link: [github.com/mahdims/2D-Bin-packing](https://github.com/mahdims/2D-Bin-packing)


