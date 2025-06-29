// compare two frequency tables to see which frequencies are present in one and not the other

const fs = require("fs")

const FILE_PATH_1 = "./freq-results/cv_cru_37_38.csv"
const FILE_PATH_2 = "./freq-results/cv_cru_64_65.csv"

const splitRe = /,(?=(?:(?:[^"]*"){2})*[^"]*$)/ // regex to match commas not in parentheses

const data1 = fs.readFileSync(FILE_PATH_1, "utf8").split("\n").map(l => l.split(splitRe).map(v => isNaN(v) ? v : Number(v)))
const data2 = fs.readFileSync(FILE_PATH_2, "utf8").split("\n").map(l => l.split(splitRe).map(v => isNaN(v) ? v : Number(v)))

const freqCol = data1[0].indexOf("Frequency")

data1.shift()
data2.shift()

console.log(data1)

// match the frequencies from one file to another
for (var i = 0; i < data1.length; i++) {
	for (var ii = 0; ii < data2.length; ii++) {
		const f1 = data1[i][freqCol], f2 = data2[ii]?.[freqCol]
		if (f2 == undefined) continue
		if (Math.abs(f1-f2) <= 0.08 && pctDiff(f1, f2) <= 0.33) {
			// if two frequencies are close enough we can "match them together"
			data1[i] = undefined
			data2[ii] = undefined
			break
		}
	}
}


// print the unmatched frequencies now
console.log("Frequencies in table 1 but not table 2")
console.log(data1.map((v, i) => v === undefined ? v : [i+1, v[freqCol]]).filter(v => v !== undefined))

console.log("\n\nFrequencies in table 2 but not table 1")
console.log(data2.map((v, i) => v === undefined ? v : [i+1, v[freqCol]]).filter(v => v !== undefined))



function pctDiff(x, y) {
	return 2*Math.abs(x-y)/(x+y)
}
