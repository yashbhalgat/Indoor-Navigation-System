myfile = open('test_cases_for_Quest.csv','r')
s = myfile.read()
myfile.close()
print s
f = open('sample.txt','w')
f.write(s)
f.close()
