import glob
import csv
import numpy as np

with open("Datas/diff_of_diffs.csv", 'w') as f:
    writer = csv.writer(f)
    file_list = glob.glob('Datas/1e-*.csv')
    for file_name in file_list:
        print(file_name)
        with open(file_name) as f:
            reader = csv.reader(f)
            datas = [row for row in reader]

        new_datas = []
        for data in datas:
            data = [float(i) for i in data]
            new_datas.append(data)
        new_datas = np.array(new_datas)
        result = np.zeros(6)
        for data in new_datas:
            result += data[12:]

        result /= len(new_datas)
        print(result)
        writer.writerow([file_name, *result])
