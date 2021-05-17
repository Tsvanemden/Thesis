#To run this program you'll need the matplotlib library installed? and a textfile containg the genomic sequence of a chromosome as a text (not fasta format only atgc's etc).
#A wee note: this script will show you a plot, x out of the plot to continue to work in you command line
import re
import matplotlib.pyplot as plt
chr = open("pNSU70.txt")
chr_contents = chr.read()
chr_contents = str(chr_contents)
chr_contents = chr_contents.replace("\n","")
chr_contents = chr_contents.upper()

print(len(chr_contents))

chr_contents_PolyAT = re.sub(r"[AT]{20}", "XXXXXXXXXXXXXXXXXXXX", chr_contents)
chr_contents_PolyAT = re.sub(r"[AT]{19}", "XXXXXXXXXXXXXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{18}", "XXXXXXXXXXXXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{17}", "XXXXXXXXXXXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{16}", "XXXXXXXXXXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{15}", "XXXXXXXXXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{14}", "XXXXXXXXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{13}", "XXXXXXXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{12}", "XXXXXXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{11}", "XXXXXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{10}", "XXXXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{9}", "XXXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{8}", "XXXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{7}", "XXXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{6}", "XXXXXX", chr_contents_PolyAT)
chr_contents_PolyAT = re.sub(r"[AT]{5}", "XXXXX", chr_contents_PolyAT)

#You wanna check every ... bases what the AT content is
steps = 1

#The ... many bp up and down stream, basically the size of the region you wanna check/ 2
span = 50

#From where to start. start from the variable steps, or any number you like
start  = 50
start_PolyAT = 50

#your stop can be either the end of the chromosome (the len() thing) or any bp.
#stop = len(chr_contents)
stop = 7226

#create empty lists that will hold the x and y coordinates of our graph
x = []
y = []
z = []

#This loop, for as long a your start isn't past your stop, first gets your ROI (region o' intrest)
#Then calculates the AT content and adds that to the lists.
#Then changes start to move to the next point and break the loop eventually
while start < stop:
    ROI = chr_contents[(start-span):(start+span)]
    AT_content = (ROI.count("A") + ROI.count("T"))/len(ROI)*100  
    x.append(start)
    y.append(AT_content)
    start = start + steps

while start_PolyAT < stop:
    ROI = chr_contents_PolyAT[(start_PolyAT-span):(start_PolyAT+span)]
    AT_content_PolyAT = ((ROI.count("X"))/len(ROI)*100)  
    z.append(AT_content_PolyAT)
    start_PolyAT = start_PolyAT + steps

#print (chr_contents[2750:3250])
#plot the s and x (bp number) and y (AT_content) and the number of PolyA streches 
#ax.spines["top"].set_visible(False)
#ax.spines["right"].set_visible(False)
plt.figure(1)
plt.subplot(211)
plt.ylim( 0, 100)
plt.xlim(150, stop)
plt.plot(x,y, color = "k", label = "A/T Content", linewidth = 1)
#plt.plot(x,z, label = "% A/T Tracts")
#plt.axvspan(2964, 3069, color='g', alpha=0.5, lw=0,label = "dh")
#plt.axvspan(2339, 2449, color='b', alpha=0.5, lw=0,label = "dg")
#plt.axvspan(30813, 36448, color='y', alpha=0.5, lw=0,label = "imr")
#plt.axvspan(36449, 40555, color='r', alpha=0.5, lw=0,label = "cnt1")
plt.axvspan(40567, 46202, color='y', alpha=0.5, lw=0)
plt.axvspan(46550, 50480, color='b', alpha=0.5, lw=0)
plt.axvspan(51563, 56352, color='g', alpha=0.5, lw=0)
plt.axhline(y=63.94, color='#666666', linestyle='--')
plt.legend()
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.suptitle('Sequence analysis of Telomere', fontsize=12)
plt.ylabel("A/T Content in %")
#plt.xlabel("Distance from telomere")

plt.subplot(212)
plt.ylim( 0, 100)
plt.xlim(0, stop)
#plt.plot(x,y, label = "% A/T Content")
plt.plot(x,z, color = "k", label = "A/T Tracts", linewidth = 1)
#plt.axvspan(21498, 25453, color='g', alpha=0.5, lw=0,label = "dh")
#plt.axvspan(26536, 30466, color='b', alpha=0.5, lw=0,label = "dg")
#plt.axvspan(30813, 36448, color='y', alpha=0.5, lw=0,label = "imr")
#plt.axvspan(36449, 40555, color='r', alpha=0.5, lw=0,label = "cnt1")
plt.axvspan(40567, 46202, color='y', alpha=0.5, lw=0)
plt.axvspan(46550, 50480, color='b', alpha=0.5, lw=0)
plt.axvspan(51563, 56352, color='g', alpha=0.5, lw=0)
plt.legend()
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.suptitle('Sequence analysis of Telomere', fontsize=12)
#plt.xticks([0, 5000], [0, "", 2.5, "", 5])
plt.ylabel("A/T Tracts in %")
plt.xlabel("Distance from telomere in bases")
plt.savefig("at_cont_pNSU70.ps")

plt.show()

# s.append(start)
#   t.append(PolyAT)
#PolyAT = len(re.findall(r"A{5,100}", ROI)) + len(re.findall(r"T{5,100}", ROI))

#plt.axvspan(41, 789, color='y', alpha=0.5, lw=0, label = "TAS")
#plt.axvspan(2086, 5532, color='y', alpha=0.5, lw=0,)
#plt.axvspan(6274, 6966, color='y', alpha=0.5, lw=0, )
#plt.axvspan(10088, 15749, color='y', alpha=0.5, lw=0,)
#plt.axvspan(17706, 19361, color='y', alpha=0.5, lw=0, )



