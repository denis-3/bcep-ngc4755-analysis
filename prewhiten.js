// analysis tool combining lomb scargle and least squares prewhitening

const fs = require("fs")
const { plot } = require("nodeplotlib");

// for the freq range the 1/T limit is 1/52.4905 for 2 secotrs
// oversample by 5 means

const DATA_FILE_NAME = "./tess-data/bw_cru_37_38.csv"
const PW_ITRS = 9999 // prewhiten itrs
const FREQ_RANGE = [0.1, 10] // frequency range
const LOMB_TRIES = 10_394 // Lomb-Scargle tries [this oversamples by 20]
const SQUARES_TRIES = 10_000 // how many points to try for phase
const NOISE_INT_WIDTH = 0.075 // width of the internal in which to find noise [days^-1]

//
// // Get The data
//

const parsedFiled = parseDataFile(DATA_FILE_NAME)
const data = parsedFiled.data
const DATA_Y_MEAN = parsedFiled.yMean

// Lomb-Scargle Least-Squares prewhiten
function lslsPrewhiten() {
	const periodogramX = []
	const periodogramYs = []
	const tableRows = [] // {freqNum, snr, freq, amp}
	var prevAmp = 0
	var prevFreq = 0
	for (var i = 0; i < PW_ITRS; i++) {
		// lomb scargle first
		var bestFreq = 0
		var bestPower = 0
		const periodogramY = []
		var nextLogTime = Date.now() + 2000
		for (var ii = 1; ii < LOMB_TRIES; ii++) {
			const thisFreq = FREQ_RANGE[0] + (FREQ_RANGE[1] - FREQ_RANGE[0]) * ii / LOMB_TRIES

			const thisTau = lombScargleTau(thisFreq, data[0])
			const thisCosPart = lombScargleCosPart(thisFreq, data[0], data[1], thisTau)
			const thisSinPart = lombScargleSinPart(thisFreq, data[0], data[1], thisTau)
			const thisLombScarglePower = (thisCosPart + thisSinPart) / 2

			if (thisLombScarglePower > bestPower) {
				bestFreq = thisFreq
				bestPower = thisLombScarglePower
			}

			if (i == 0) periodogramX.push(thisFreq)
			periodogramY.push(2 * Math.sqrt(thisLombScarglePower / data[0].length))

			if (Date.now() >= nextLogTime) {
				console.log("Percent done with Lomb-Scargle", Math.round(100 * ii / LOMB_TRIES))
				nextLogTime += 2000
			}
		}

		periodogramYs.push(periodogramY)

		// calculate SNR for the previous freq
		if (i > 0) {
			// code to do noise interal
			// interval start i of the noise interval
			// const intStartI = prevFreq > FREQ_RANGE[0] + NOISE_INT_WIDTH/2
			// 	? (prevFreq - NOISE_INT_WIDTH/2 - FREQ_RANGE[0]) * LOMB_TRIES / (FREQ_RANGE[1] - FREQ_RANGE[0])
			// 	: 0
			// const intEndI = prevFreq < FREQ_RANGE[1] - NOISE_INT_WIDTH/2
			// 	? (prevFreq + NOISE_INT_WIDTH/2 - FREQ_RANGE[0]) * LOMB_TRIES / (FREQ_RANGE[1] - FREQ_RANGE[0])
			// 	: data[0].length
			// const prevSnr = prevAmp / arrayStdev(periodogramY.slice(Math.floor(intStartI), Math.ceil(intEndI+1)), 0)
			const prevSnr = prevAmp / arrayStdev(periodogramY, 0)
			console.log("signal to noise ratio:", prevSnr)
			if (prevSnr < 5.124) {
				tableRows.pop()
				break
			} else {
				tableRows[tableRows.length - 1].snr = roundToDecimal(prevSnr, 2)
			}
		}

		if (i >= 3) periodogramYs.pop()

		const amp = 2 * Math.sqrt(bestPower / data[0].length)
		console.log("best current freuqny", bestFreq, "power", bestPower)

		var bestPhase = 0
		var bestRss = -1
		// now least squares for phase
		for (var iii = 0; iii < SQUARES_TRIES; iii++) {
			const thisPct = iii / SQUARES_TRIES
			const thisPhase = 2*Math.PI * thisPct
			const thisRss = residualsSumOfSquares([amp, 2*Math.PI*bestFreq, thisPhase], data[0], data[1])
			if (thisRss < bestRss || bestRss == -1) {
				bestRss = thisRss
				bestPhase = thisPhase
			}
		}
		const bestSineParams = [amp, 2*Math.PI*bestFreq, bestPhase]
		console.log("final params, amp, freq, phase")
		console.log(amp)
		console.log(bestFreq)
		console.log(bestPhase)

		const thisRmsd = Math.sqrt(bestRss / data[0].length)
		const ampUncrt = Math.sqrt(2 / data[0].length) * thisRmsd
		const phaseUncrt = ampUncrt / bestSineParams[0]
		const freqUncrt = phaseUncrt * Math.sqrt(3) / (Math.PI * data[0][data[0].length-1])

		console.log("freq uncertainity", freqUncrt)

		data[1] = data[1].map((y, ii) => y - getSinePrediction(bestSineParams, data[0][ii]))

		tableRows.push({
			freqNum: i + 1,
			freq: roundToDecimal(bestFreq, 5),
			amp: `${roundToDecimal(amp, 1)} ± ${roundToDecimal(ampUncrt, 1)}`,
			mmagAmp: roundToDecimal(1250 * Math.log10((DATA_Y_MEAN + amp)/(DATA_Y_MEAN - amp)), 2)
		})

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
		prevAmp = amp
		prevFreq = bestFreq
	}
	console.log("total of", tableRows.length, "frequencies")
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
		// config: {staticPlot: true}
	})

	console.table(tableRows, ["freqNum", "freq", "snr", "amp", "mmagAmp"])
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

// the first fraction in lomb scargle
function lombScargleCosPart(freq, timeStamps, fluxValues, tau) {
	var numerator = 0
	var denominator = 0
	for (var i = 0; i < timeStamps.length; i++) {
		const cosArgument = 2 * Math.PI * freq * (timeStamps[i] - tau)
		numerator += fluxValues[i] * Math.cos(cosArgument)
		denominator += (Math.cos(cosArgument) ** 2)
	}
	numerator = numerator ** 2
	return numerator / denominator
}

// the second fraction in lomb scargle
function lombScargleSinPart(freq, timeStamps, fluxValues, tau) {
	var numerator = 0
	var denominator = 0
	for (var i = 0; i < timeStamps.length; i++) {
		const sinArgument = 2 * Math.PI * freq * (timeStamps[i] - tau)
		numerator += fluxValues[i] * Math.sin(sinArgument)
		denominator += (Math.sin(sinArgument) ** 2)
	}
	numerator = numerator ** 2
	return numerator / denominator
}

// tau function in lomb scargle
function lombScargleTau(freq, timeStamps) {
	var sinSum = 0
	var cosSum = 0
	for (var i = 0; i < timeStamps.length; i++) {
		const sinCosArgument = 4 * Math.PI * freq * timeStamps[i]
		sinSum += Math.sin(sinCosArgument)
		cosSum += Math.cos(sinCosArgument)
	}
	const arctanResult = Math.atan2(sinSum, cosSum)
	return arctanResult / (4 * Math.PI * freq)
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
