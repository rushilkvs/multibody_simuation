import re
s = input("File name with extension : ")
f = open(s, 'r')
f1 = open("inp.txt", 'w')
s = f.read()
pattern = re.compile(r'(= )(\S+)')
for a,b in re.findall(pattern, s):
    print(b)
    f1.write(b + '\n')

pattern = re.compile(r'\d+\.\d+\s\d+\.\d+\s\d+\.\d+')
for t in re.findall(pattern, s):
    f1.write(t + '\n')