#!/usr/bin/env python

import sys

def table(datalist):
	filenames = [i[1] for i in datalist]
	print("<tr>")
	print("<td>Filename</td>")
	for element in filenames:
		print(f"<td>{element}</td>")
	print("</tr>")

	seq_len = [i[2] for i in datalist]
	print("<tr>")
	print("<td>Sequence length</td>")
	for element in seq_len:
		print(f"<td>{element}</td>")
	print("</tr>")

	total_seq = [i[3] for i in datalist]
	print("<tr>")
	print("<td>Total Sequences</td>")
	for element in total_seq:
		print(f"<td>{element}</td>")
	print("</tr>")

	gc_content = [i[4] for i in datalist]
	print("<tr>")
	print("<td>%GC</td>")
	for element in gc_content:
		print(f"<td>{element}</td>")
	print("</tr>")
		
	return

file = sys.argv[1]

pre_list = []
post_list = []


with open(file,"r") as infile:
	infile = infile.readlines()
	infile = [element.replace("\n","") for element in infile]
	infile = [element.split(",") for element in infile]
	infile.sort(key=lambda x:(x[0],x[1]))


# 0: prepost, 1: filename, 2: seqlen, 3: nseqs, 4: gc_content, 5: path
for sample in infile:
	samplename = sample[0]
	end = sample[1]
	prepost = sample[2]
	filename = sample[3]
	seqlen = sample[4]
	nseqs = sample[5]
	gc_content = sample[6]
	path = sample[7]

	if prepost == "pre":
		pre_list.append([samplename,filename,seqlen,nseqs,gc_content,path])

	if prepost == "post":
			post_list.append([samplename,filename,seqlen,nseqs,gc_content,path])

print("<table class='table'>\
		 			<tr>\
            <thead>\
              <th>Sample</th>")

for data in pre_list:
	print(f"<th>{data[0]}</th>")

print("</thead><tbody></tr>")
print(f"<tr><td colspan='{len(pre_list)+1}' class='info'>Pre-Filter</td></tr>")
table(pre_list)

print(f"</tr><tr><td colspan='{len(post_list)+1}'class='info'>Post-Filter</td></tr>")
table(post_list)

print("</tbody></table></div><div><br>")

print ("<p>Pre-Filter Reports:</p><ul>")
for data in pre_list:
	print(f"<li><a target='_blank' href='{data[-1]}'>{data[0]}:{data[1]}</a></li>")
print("</ul>")

print("<p>Post-Filter Reports:</p><ul>")
for data in post_list:
	print(f"<li><a target='_blank' href='{data[-1]}'>{data[0]}{data[1]}</a></li>")
print("</ul>")






