import ast

readout = ast.literal_eval(open('../hb_presence_2.txt', 'r').readlines()[0].strip('\n'))
print(len(readout['2']))
