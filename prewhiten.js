// analysis tool combining lomb scargle and least squares prewhitening

const fs = require("fs")
const { plot } = require("nodeplotlib")
const { Worker } = require("worker_threads")

const DATA_FILE_NAME = "./tess-data/bw_cru_37_38.csv"
const PW_ITRS = 9999 // prewhiten it
const FREQ_RANGE = [0.1, 10] // frequency range
const LOMB_TRIES = 10_394 // Lomb-Scargle tries [this oversamples by 20]
const SQUARES_TRIES = 10_000 // how many points to try for phase
const NOISE_INT_WIDTH = 0.075 // width of the internal in which to find noise [days^-1]
const WORKER_COUNT = 10

if (WORKER_COUNT > require("os").cpus().length) {
	console.warn("\n~~~\nWarning: Using more workers than CPU cores is not advised since there will be multiple workers in a core!\n~~~\n")
}

//
// // Get The data
//

const parsedFiled = parseDataFile(DATA_FILE_NAME)
const data = parsedFiled.data
const DATA_Y_MEAN = parsedFiled.yMean
const xDataBuf = new SharedArrayBuffer(Float64Array.BYTES_PER_ELEMENT * data[0].length)
const yDataBuf = new SharedArrayBuffer(Float64Array.BYTES_PER_ELEMENT * data[0].length)
const xData = new Float64Array(xDataBuf)
const yData = new Float64Array(yDataBuf)

for (var i = 0; i < data[0].length; i++) {
	xData[i] = data[0][i]
	yData[i] = data[1][i]
}

console.log("Creating workers...")
const workers = []
for (var i = 0; i < WORKER_COUNT; i++) {
	const w = new Worker(__dirname + "/prewhiten-worker.js", { workerData: {x: xDataBuf, y: yDataBuf, freqRange: FREQ_RANGE, lombTries: LOMB_TRIES}})
	const wi = i // worker index
	w.on("error", (e) => {console.log("Worker", wi, "errored:", e)})
	w.on("exit", (s) => {console.log("Worker", wi, "exited with status code:", s)})
	workers.push(w)
}

// Lomb-Scargle Least-Squares prewhiten
async function lslsPrewhiten() {
	var periodogramX = new Array(workers.length)
	const periodogramYs = []
	const tableRows = [] // {freqNum, snr, freq, amp}
	for (var i = 0; i < PW_ITRS; i++) {
		// lomb scargle first
		var periodogramY = new Array(workers.length)
		var workersDone = 0
		console.log("Notifying workers for iteration", i, "...")
		for (var ii = 0; ii < workers.length; ii++) {
			workers[ii].removeAllListeners("message")
			workers[ii].postMessage({command: "start periodogram",
				startI: Math.floor(ii*LOMB_TRIES/workers.length),
				endI: Math.ceil((ii+1)*LOMB_TRIES/workers.length),
				iteration: i})
			const wi = ii // worker index
			workers[ii].on("message", (m) => {
				if (m.message != "finished periodogram")
				console.log()
				if (m.x.length != 0) {
					periodogramX[wi] = m.x
				}
				periodogramY[wi] = m.y
				workersDone++
			})
		}

		var nextLogTime = Date.now() + 2000
		while (workersDone < workers.length) {
			await millis(250)
			if (Date.now() >= nextLogTime) {
				console.log("Workers done:", workersDone, "out of", workers.length)
				nextLogTime += 2000
			}
		}
		console.log("Workers finished!")

		if (i == 0) periodogramX = periodogramX.flat()
		periodogramY = periodogramY.flat()

		var bestFreq = 0
		var bestAmp = 0
		for (var ii = 0; ii < periodogramX.length; ii++) {
			if (periodogramY[ii] > bestAmp) {
				bestFreq = periodogramX[ii]
				bestAmp = periodogramY[ii]
			}
		}

		// significance check
		const thisSnr = bestAmp / arrayMedian(periodogramY)
		console.log("Signal to Noise Ratio:", thisSnr)
		if (thisSnr <= 5.124) {
			periodogramYs.push(periodogramY) // residual amplitude spectrum
			break
		}

		if (i < 3) periodogramYs.push(periodogramY)

		console.log("Finding phase...")
		var bestPhase = 0
		var bestRss = -1
		// now least squares for phase
		for (var iii = 0; iii < SQUARES_TRIES; iii++) {
			const thisPhase = 2*Math.PI * iii / SQUARES_TRIES
			const thisRss = residualsSumOfSquares([bestAmp, 2*Math.PI*bestFreq, thisPhase], xData, yData)
			if (thisRss < bestRss || bestRss == -1) {
				bestRss = thisRss
				bestPhase = thisPhase
			}
		}
		const bestSineParams = [bestAmp, 2*Math.PI*bestFreq, bestPhase]
		console.log("Final amplitude, frequency, and phase")
		console.log(bestAmp)
		console.log(bestFreq)
		console.log(bestPhase)

		const thisRmsd = Math.sqrt(bestRss / xData.length)
		const ampUncrt = Math.sqrt(2 / xData.length) * thisRmsd
		const phaseUncrt = ampUncrt / bestSineParams[0]
		const freqUncrt = phaseUncrt * Math.sqrt(3) / (Math.PI * (xData[xData.length-1] - xData[0]))

		for (var ii = 0; ii < yData.length; ii++) {
			yData[ii] = yData[ii] - getSinePrediction(bestSineParams, xData[ii])
		}

		tableRows.push({
			freqNum: i + 1,
			freqDisplay: `${roundToDecimal(bestFreq, 5)} ± ${roundToDecimal(freqUncrt, 5)}`,
			snr: roundToDecimal(thisSnr, 2),
			ampDisplay: `${roundToDecimal(bestAmp, 1)} ± ${roundToDecimal(ampUncrt, 1)}`,
			mmagAmp: roundToDecimal(1250 * Math.log10((DATA_Y_MEAN + bestAmp)/(DATA_Y_MEAN - bestAmp)), 2),
			freq: roundToDecimal(bestFreq, 5),
			amp: roundToDecimal(bestAmp, 1),
		})

		// for debug
		if (i == 10000) {
			const trace = {
				x: periodogramX,
				y: periodogramY,
				mode: "markers",
				marker : {
					color: "green",
					size: 4
				},
				name: "Periodogram after iteration " + String(i+1)
			}

			const layout = {
				title: "periodogram",
				width: 700,
				height: 375
			}

			plot({
				data: [trace],
				layout
			})
			break
		}
		console.log("\n")
	}
	console.log("Shutting down workers...")
	workers.forEach(w => w.postMessage({command: "exit"}))


	console.log("Total of", tableRows.length, "frequencies")
	const layout = {
		showlegend: false,
		xaxis: {
			tickfont: {
				size: 20,
			},
			dtick: 1,
			tick0: 0
		},
		yaxis: {
			tickfont: {
				size: 20
			},
			tick0: 0
		},
		grid: {
			rows: 2,
			columns: 2,
			pattern: "independent"
		},
		width: 700*2,
		height: 300*2,
		annotations: [],
		margin: { t: 50 }
	}
	const traces = []
	for (var i = 0; i < 4; i++) {
		traces.push({
			x: periodogramX,
			y: periodogramYs[i],
			mode: "lines",
			line: {
				color: "black",
				width: 1.66,
				shape: "spline"
			},
			xaxis: `x${i+1}`,
			yaxis: `y${i+1}`,
		})
		layout[`xaxis${i+1}`] = JSON.parse(JSON.stringify(layout.xaxis))
		layout[`yaxis${i+1}`] = JSON.parse(JSON.stringify(layout.yaxis))

		layout.annotations.push({
			xref: `x${i+1} domain`,
			yref: `y${i+1} domain`,
			x: 0.08,
			y: 0.8,
			text: i < 3 ? `<i>f<sub>${i+1}</sub></i>=${roundToDecimal(tableRows[i].freq, 3)} d<sup>-1</sup>` : "Residuals",
			xanchor: "left",
			yanchor: "bottom",
			showarrow: false,
			font: {
				size: 22
			},
			bgcolor: "white"
		})
	}

	layout.yaxis.title = {text: 'Amplitude (counts)', standoff: 40, font: {size:20}}
	layout.yaxis.automargin = true

	layout.yaxis3.title = {text: 'Amplitude (counts)', standoff: 40, font: {size:20}}
	layout.yaxis3.automargin = true

	layout.xaxis3.title = 'Frequency (days<sup>-1</sup>)'
	layout.xaxis3.titlefont = {size: 20}

	layout.xaxis4.title = 'Frequency (days<sup>-1</sup>)'
	layout.xaxis4.titlefont = {size: 20}

	plot({
		data: traces,
		layout,
	})

	console.table(tableRows, ["freqNum", "freqDisplay", "snr", "ampDisplay", "mmagAmp"])

	const csvSave = "./freq-results/" + DATA_FILE_NAME.split("/").reverse()[0]
	fs.writeFileSync(csvSave,
		"Frequency,SNR,Amplitude,MmagAmplitude" +
		tableRows.map(r =>
			`\n"f${r.freqNum}, ${r.freqDisplay}",${r.snr},${r.ampDisplay},${r.mmagAmp}`
		).join(""))
	console.log("Saved frequencies to", csvSave)
}

lslsPrewhiten()


//
// // Helper funcs
//

function arrayMean(arr) {
	if (arr.length == 0) throw Error("Zero-length array")
	return arr.reduce((a, c) => a + c, 0) / arr.length
}

function arrayStdev(arr, arrMean) {
	if (arr.length == 0) throw Error("Zero-length array")
	const squareSum = arr.reduce((a, c) => a + (c - arrMean) ** 2, 0)
	return Math.sqrt(squareSum / arr.length)
}

function arrayMedian(values) {
	const sorted = values
		.slice()
		.sort((a, b) => a - b);

	const len = sorted.length;
	const mid = Math.floor(len / 2);

	if (len % 2 == 0) {
		return (sorted[mid - 1] + sorted[mid]) / 2;
	} else {
		return sorted[mid];
	}
}

function parseDataFile(fileName) {
	const initData = fs.readFileSync(fileName, "utf8").split("\n")
	.map(row =>
		row.split(",").map(
			num => Number(num)
		)
	)

	initData.shift()

	const thisData = [
		[],
		[]
	]

	for (var i = 0; i < initData.length; i++) {
		if (isNaN(initData[i][0]) || isNaN(initData[i][1])) continue;

		thisData[0].push(initData[i][0])
		thisData[1].push(initData[i][1])
	}

	const dataYMean = arrayMean(thisData[1])
	thisData[1] = thisData[1].map(y => y - dataYMean) // subtract mean from data

	return {data: thisData, yMean: dataYMean}
}

function getSinePrediction(sineParams, x) {
	return sineParams[0] * Math.sin(sineParams[1] * x + sineParams[2])
}

function residualsSumOfSquares(sineParams, inputX, inputY) {
	var squareSum = 0
	// find sum of squares
	for (var i = 0; i < inputY.length; i++) {
		const thisPrediction = getSinePrediction(sineParams, inputX[i])
		squareSum += (inputY[i] - thisPrediction) ** 2
	}

	return squareSum
}

function roundToDecimal(inp, dec) {
	return Math.round(inp * (10**dec)) / (10**dec)
}

async function millis(m) {
	return new Promise((resolve, reject) => {
		setTimeout(() => {
			resolve(true)
		}, m)
	})
}
