'''
    B4TM patient - label variables orderer from file
    Created on 04-04-2018
    Author: Mark Grannetia
    mga 202, 2624692
'''

# This program distracts the variables patients and their associated label (subgroup).
SAMPLE_AND_DATA_FILE = 'Train_call.txt'


def process_file(input_file):
    for line in input_file:
        # combination = line.strip('\n').replace('\"','').split('\t')
        # # Cleans-up txt and splits combinations of patient and subgroup
        # sample = combination[0]
        # subgroup = combination[-1]
        # print('%s \t %s' % (sample, subgroup))
        print(line)


samples_and_data = open(SAMPLE_AND_DATA_FILE).readlines()

print(samples_and_data)

process_file(samples_and_data)


# print 'The labels are: '
# for student in students:
#     print process_names(student)
# print 'End of report'














































































