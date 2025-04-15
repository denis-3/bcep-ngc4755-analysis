// An analysis tool that can do least-squares prewhitening fits, Lomb-Scargle periodogram, and O-C graphs
const fs = require("fs")
const crypto = require("crypto")
const { plot } = require("nodeplotlib");
const prando = require("prando")
const regression = require("regression")

// General config
const ANALYSIS_TYPE = "prewhiten" // "prewhiten," "period," or "oc"
const MAX_DATA_POINTS = 10e6 // Max data points to use (more means laggy, but more accurate results)
const STAR_NAME = "" // name of star (BW Cru or CV Cru)
const DATA_FILE_NAME = "" // data file name
const RAND_SEED = md5Hash(String(Math.random())) // init random seed
const GLOBAL_RNG = new prando(RAND_SEED)


// Prewhitening config
const PREWHITEN_RAND_GUESS = 100_000 // amount of rand guesses before main prewhiten
const PREWHITEN_CT = 0 // [minimum] Prewhiten depth
const PREWHITEN_MAX_CT = 10
const PREWHITEN_TRIES = 400_000 // how many times each sine parameter is adjusted per depth
const PARAM_TUNE_ITRS = 7 // How many iterations of tuning
const DATA_EPOCH = 0 // beginning of data

const PREWHITEN_PERIOD_RANGE = [0.11111, 60] // range of periods fo try for prewhitn
const PREWHITEN_FREQ_RANGE = PREWHITEN_PERIOD_RANGE.map(p => 1 / p).reverse()
const PW_SNR_THRSH = 3 // SNR of prewhiten threshold
const PW_BIN_SIZE = 10 // bin size for noise calc


// Period-folding config
const PERIOD_RANGE = [0.07, 60] // range of periods to try for period folding
const PERIOD_INCTS = 30000 // increments in the frequencies that are tried

PERIOD_RANGE[0] = PERIOD_RANGE[0] * 60 * 60 * 24
PERIOD_RANGE[1] = PERIOD_RANGE[1] * 60 * 60 * 24

const FREQ_RANGE = [1 / PERIOD_RANGE[1], 1 / PERIOD_RANGE[0]]


// O-C config (OC can support two stars; the config is per star usually)
const OTHER_FILE_NAME = ""
const MAIN_FREQUENCIES = [1, 2] // the frequency for making chunks (days^-1) for both stars
const SINE_FIT_GUESSES = 50_000 // count of initial random guesses
const OC_TUNE_TRIES = 100_000 // how many times each param gets tuned in total
const OC_TUNE_ITRS = 4 // how many total tuning iterations

//
// // Get The data
//

const data = parseDataFile(DATA_FILE_NAME)

const DATA_Y_MEAN = arrayMean(data[1])
data[1] = data[1].map(y => y - DATA_Y_MEAN) // subtract mean from data

if (ANALYSIS_TYPE == "prewhiten") {
	console.log("Using seed:", RAND_SEED)
	prewhitenAnalysis()
} else if (ANALYSIS_TYPE == "period") {
	data[0] = data[0].map(x => x * 60 * 60 * 24) // convert MJD to seconds
	periodAnalysis()
} else if (ANALYSIS_TYPE == "oc") {
	ocGraphs()
} else {
	console.log("Unknown analysis", ANALYSIS_TYPE)
}

//
// // Period folding and periodogram analysis
//

function periodAnalysis() {
	var bestPeriod = -1
	var bestLombScarglePower = -1
	var bestPeriodFoldedX = []
	var periodogramData = [
		[],
		[]
	]
	console.log("~~~Finding best period~~~")
	const periodFoldedData = JSON.parse(JSON.stringify(data))
	var nextLogTime = Date.now() + 2000
	for (var i = 1; i < PERIOD_INCTS; i++) {
		const thisFreq = FREQ_RANGE[0] + (FREQ_RANGE[1] - FREQ_RANGE[0]) * i / PERIOD_INCTS

		const thisTau = lombScargleTau(thisFreq, periodFoldedData[0])
		const thisCosPart = lombScargleCosPart(thisFreq, periodFoldedData[0], periodFoldedData[1], thisTau)
		const thisSinPart = lombScargleSinPart(thisFreq, periodFoldedData[0], periodFoldedData[1], thisTau)
		const thisLombScarglePower = (thisCosPart + thisSinPart) / 2

		if (thisLombScarglePower > bestLombScarglePower) {
			const thisPeriod = 1 / thisFreq
			bestLombScarglePower = thisLombScarglePower
			bestPeriod = thisPeriod
			bestPeriodFoldedX = periodFoldedData[0].map(x => (x / thisPeriod) % 1)
		}
		periodogramData[0].push(thisFreq) // push frequency as x axis to have frequency-based graph
		periodogramData[1].push(thisLombScarglePower)
		if (Date.now() >= nextLogTime) {
			console.log("Percent done", Math.round(100 * i / PERIOD_INCTS))
			nextLogTime += 2000
		}
	}

	console.log("The best period is (hours)", (bestPeriod / 3600), "frequency", 86400 / bestPeriod, "with lomb scargle", bestLombScarglePower)

	periodogramData[1] = periodogramData[1].map(y => y * 100 / bestLombScarglePower)
	data[0] = data[0].map(x => x / (60*60*24)) // map seconds back to days

	const layout = {
		showlegend: false,
		staticPlot: true,
		margin: { t: 50 },
		xaxis2: {
			title: 'Phase',
			titlefont: {
				size: 20
			},
			tickfont: {
				size: 20
			}
		},
		yaxis2: {
			title: 'Flux (counts)',
			titlefont: {
				size: 20
			},
			tickfont: {
				size: 20
			}
		},
		xaxis: {
			title: 'Frequency (days<sup>-1</sup>)',
			titlefont: {
				size: 20
			},
			tickfont: {
				size: 20
			}
		},
		yaxis: {
			title: 'Lomb-Scargle Power',
			titlefont: {
				size: 20
			},
			tickfont: {
				size: 20
			}
		},
		width: 704*2,
		height: 435,
		grid: {
			rows: 1,
			columns: 2,
			pattern: "independent"
		}
	}

	const periodFoldedTrace = {
		x: bestPeriodFoldedX,
		y: data[1],
		mode: "markers",
		marker: {
			autocolorscale: false,
			colorscale: "Portland",
			color: data[0],
			size: 5,
			opacity: 0.3,
			showscale: true,
			colorbar: {
				x: 1.04,
				title: {
					text: "Time (Days)",
					font: { size: 20 }
				},
				tickfont: {
					size: 20
				}
			}
		},
		name: "TESS Data",
		xaxis: "x2",
		yaxis: "y2",
	}

	const periodogramTrace = {
		x: periodogramData[0].map(x => x * 86400),
		y: periodogramData[1],
		mode: "lines",
		line: {
			color: "black",
			width: 3,
			shape: "spline"
		},
		name: "Lomb-Scargle Power",
		xaxis: "x1",
		yaxis: "y1"
	}

	plot({
		data: [periodogramTrace, periodFoldedTrace],
		layout,
		config: {staticPlot: true}
	});
}

//
// // Prewhiten Anlysis function
//

function prewhitenAnalysis() {
	const dataMean = arrayMean(data[1])
	const dataRange = Math.max(...data[1]) - Math.min(...data[1])
	var bestSquareSum = -1;
	// initial guess to sine fit
	const guessSineParams = [
		dataRange, // amplitude
		2 * Math.PI * PREWHITEN_FREQ_RANGE[0], // b
		Math.PI, // phase angle
		0 // vertical shift
	]
	var bestSineParams = JSON.parse(JSON.stringify(guessSineParams))
	var latestPwResids = []
	const finalizedSineParams = []
	const sinesMetaData = [] // array of arrays [SNR, mmag amplitude, phase at epoch, mmag amplitude uncertainty, phase uncert, freq uncert]
	const allPrewhitenResids = []

	const nanData = JSON.parse(JSON.stringify(data))
	if (STAR_NAME == "BW Cru") {
		for (var i = 0; i < nanData[0].length; i++) {
			if (nanData[0][i] > 12) {
				nanData[0].splice(i, 0, 12)
				nanData[1].splice(i, 0, NaN)
				break
			}
		}
	} else if (STAR_NAME == "CV Cru") {
		const nanX = [7.3, 14.8, 21.6]
		var nanI = 0
		for (var i = 0; i < nanData[0].length; i++) {
			if (nanData[0][i] > nanX[nanI]) {
				nanData[0].splice(i, 0, 12)
				nanData[1].splice(i, 0, NaN)
				nanI += 1
				if (nanI == nanX.length) break
			}
		}
	}
	var dataToUse = JSON.parse(JSON.stringify(data))

	console.log("~~~ Main Fit Find ~~~")

	// iterate through the data
	for (var i = 0; i <= PREWHITEN_MAX_CT; i++) {
		// first the random guessing
		const dataRange = Math.max(...dataToUse[1]) - Math.min(...dataToUse[1])
		for (var ii = 0; ii < PREWHITEN_RAND_GUESS; ii++) {
			const thisSineParams = getRandomSineParams(GLOBAL_RNG, PREWHITEN_FREQ_RANGE, dataRange)
			const thisSquareSum = residualsSumOfSquares(thisSineParams, dataToUse[0], dataToUse[1])

			if (thisSquareSum < bestSquareSum || bestSquareSum == -1) {
				bestSquareSum = thisSquareSum
				bestSineParams = thisSineParams
			}
		}
		console.log("Completed initial random guesses")

		// then "fine-tuning"
		for (var ii = 0; ii < PARAM_TUNE_ITRS; ii++) {
			const newTunedSine = tuneSineParams(bestSineParams, PREWHITEN_FREQ_RANGE,
				PREWHITEN_TRIES / PARAM_TUNE_ITRS, dataToUse[0], dataToUse[1])
			if (newTunedSine.bestRss < bestSquareSum) {
				bestSquareSum = newTunedSine.bestRss
				bestSineParams = newTunedSine.params
			} else break // if there are no changes now, there will be no changes in the rest of the itrs

			console.log("Completed tune iterations", ii)
		}
		// residuals on prewhiten
		const prewhitenResiduals = dataToUse[1].map((y, v) => y - getSinePrediction(bestSineParams, dataToUse[0][v]))

		// noise first
		const thisMidLine = findMiddleLine(PW_BIN_SIZE, dataToUse[0], prewhitenResiduals)
		const thisNoise = findNoise(thisMidLine[0], thisMidLine[1], dataToUse[0], prewhitenResiduals)
		const noisePower = findNoisePower(thisNoise)

		const thisSnr = bestSineParams[0]**2 / 2 / noisePower
		const thisMmagAmp = 1250 * Math.log10((DATA_Y_MEAN + bestSineParams[0])/(DATA_Y_MEAN - bestSineParams[0]))
		const thisRmsd = Math.sqrt(bestSquareSum / dataToUse[0].length)
		const thisAmpUncrt = Math.sqrt(2 / dataToUse[0].length) * thisRmsd
		const thisPhaseUncrt = thisAmpUncrt / bestSineParams[0]
		const thisFreqUncrt = thisPhaseUncrt * Math.sqrt(3) / (Math.PI * dataToUse[0][dataToUse[0].length-1])
		const phaseAtEpoch = findPhaseAtJd(bestSineParams[1], bestSineParams[2], DATA_EPOCH)
		console.log("the snr", thisSnr)
		console.log("mmag amp", thisMmagAmp)
		console.log("phase at epoch", phaseAtEpoch)

		// analysis is done
		if (thisSnr < PW_SNR_THRSH && i > PREWHITEN_CT) {
			console.log("Insignificant sine curve on iteration", i, "\n")
			break
		}

		allPrewhitenResids.push(prewhitenResiduals)
		finalizedSineParams.push(JSON.parse(JSON.stringify(bestSineParams)))
		sinesMetaData.push([thisSnr, thisMmagAmp, phaseAtEpoch, thisAmpUncrt, thisPhaseUncrt, thisFreqUncrt])
		console.log("Best params for prewhiten", i, ": ", bestSineParams)


		bestSquareSum = -1
		dataToUse = [data[0], JSON.parse(JSON.stringify(prewhitenResiduals))]
		console.log("Completed prewhiten iteration", i, "\n")
	}

	console.log("All best sine params", finalizedSineParams)
	console.log("Receipt for star", STAR_NAME, "(seed is", RAND_SEED, ")")
	const tabledData = []
	for (var i = 0; i < finalizedSineParams.length; i++) {
		tabledData.push({
			freq: finalizedSineParams[i][1] / (2 * Math.PI),
			freqUncrt: sinesMetaData[i][5],
			snr: sinesMetaData[i][0],
			ampTessFlux: finalizedSineParams[i][0],
			ampTessFluxUncrt: sinesMetaData[i][3],
			ampMmag: sinesMetaData[i][1],
			epochPhase: sinesMetaData[i][2],
			phaseUncrt: sinesMetaData[i][4]
		})
	}
	console.table(tabledData)

	//
	// // Plot the data
	//

	console.log("~~~ Display Data ~~~")

	const rowCount = Math.ceil(finalizedSineParams.length / 2)
	const layout = {
		showlegend: false,
		margin: { t: 50 },
		xaxis: {
			title: 'Time (Days)',
			titlefont: {
				size: 20
			},
			tickfont: {
				size: 20
			},
			tick0: 0,
			dtick: 1
		},
		yaxis: {
			title: 'Flux (Counts)',
			titlefont: {
				size: 20
			},
			tickfont: {
				size: 20
			},
			tick0: 0,
			dtick: 200
		},
		width: 704*2,
		height: 435 * rowCount,
		grid: {
			rows: rowCount,
			columns: 2,
			pattern: 'independent'
		},
		annotations: []
	}

	const allTraces = []

	for (var i = -1; i < finalizedSineParams.length-1; i++) {
		// tess data or prewhitenening residuals
		// first layout config
		if (i > -1) {
			layout[`xaxis${i+2}`] = JSON.parse(JSON.stringify(layout.xaxis))
			layout[`yaxis${i+2}`] = JSON.parse(JSON.stringify(layout.yaxis))
		}
		// show the proportion of the data equal to the first six days
		var dataPct = 1
		if (Math.PI / finalizedSineParams[i+1][1] < 0.5) {
			dataPct = data[0].findIndex(t => t >= 6) / data[0].length
		} else if (i > -1) {
			delete layout[`xaxis${i+2}`].dtick
		}
		const cutoffIdx = Math.floor(data[0].length * dataPct)
		const thisYData = i == -1 ? data[1].slice(0, cutoffIdx) : allPrewhitenResids[i].slice(0, cutoffIdx)
		const thisXData = data[0].slice(0, cutoffIdx)
		const thisDataTrace = {
			x: thisXData,
			y: thisYData,
			mode: "markers",
			marker: {
				color: "grey",
				size: 2.5,
				opacity: 0.45
			},
			name: "Data" + String(i),
			xaxis: "x" + String(i+2),
			yaxis: "y" + String(i+2)
		}

		allTraces.push(thisDataTrace)

		if (i < finalizedSineParams.length-1) {
			// sine curve data - a point every 17 minutes
			const constIntervalX = [data[0][0]]
			const cIntvP = 17 // minutes (constant interval period)
			const totalConstIntervals = (data[0][data[0].length * dataPct - 1] - data[0][0]) / (cIntvP / (24 * 60))
			for (var ii = 1; ii < totalConstIntervals; ii++) {
				constIntervalX.push(ii * cIntvP / (24 * 60))
			}

			const thisSineTrace = {
				x: constIntervalX,
				y: constIntervalX.map(x => getSinePrediction(finalizedSineParams[i+1], x)),
				mode: "lines",
				line: {
					color: "rgba(0, 0, 0, 0.7)",
					width: 1.5,
					shape: "spline"
				},
				name: "Sine",
				xaxis: "x" + String(i+2),
				yaxis: "y" + String(i+2)
			}
			allTraces.push(thisSineTrace)
		}
	}

	for (var i = 0; i < finalizedSineParams.length; i++) {
		layout.annotations.push({
			xref: `x${i+1} domain`,
			yref: `y${i+1} domain`,
			x: 0.08,
			y: 0.87,
			text: `<i>f<sub>${i+1}</sub></i>`,
			xanchor: "left",
			yanchor: "bottom",
			showarrow: false,
			font: {
				size: 26
			},
			bgcolor: "white"
		})
	}

	plot({
		data: allTraces,
		layout,
		config: {staticPlot: true}
	})
}

//
// // Period change analysis
//

function ocGraphs() {
	const ocTraces = []
	const bestFitTraces = []
	var dataToUse = JSON.parse(JSON.stringify(data))
	// s for star
	for (var s = 0; s < 2; s++) {
		if (s == 0) console.log("Starting ", DATA_FILE_NAME)
		else {
			console.log("Starting ", OTHER_FILE_NAME)
			dataToUse = parseDataFile(OTHER_FILE_NAME)
		}
		var nextPeriodCutoff = 1 / MAIN_FREQUENCIES[s]
		var chunkedXData = []
		var chunkedYData = []
		const maxima = []
		var nextUpdateTime = Date.now() + 2000
		var i = 0
		while (i < dataToUse[0].length) {
			chunkedXData.push(dataToUse[0][i])
			chunkedYData.push(dataToUse[1][i])
			i++
			// do the sine fitting if we have enough points
			if (dataToUse[0][i] >= nextPeriodCutoff) {
				var bestRss = -1
				var bestSineParams = [1, 2, 3, 0]
				const chunkXStart = chunkedXData[0]
				chunkedXData = chunkedXData.map(x => x - chunkXStart)
				const thisYMean = arrayMean(chunkedYData)
				chunkedYData = chunkedYData.map(y => y - thisYMean)
				const dataRange = Math.max(...chunkedYData) - Math.min(...chunkedYData)
				// get the best sine curve from randoms
				for (var ii = 0; ii < SINE_FIT_GUESSES; ii++) {
					const thisSineParams = getRandomSineParams(GLOBAL_RNG, [MAIN_FREQUENCIES[s] * 0.5, MAIN_FREQUENCIES[s] * 2], dataRange)
					const thisSquareSum = residualsSumOfSquares(thisSineParams, chunkedXData, chunkedYData)

					if (thisSquareSum < bestRss || bestRss == -1) {
						bestRss = thisSquareSum
						bestSineParams = thisSineParams
					}
				}

				// tune it now
				for (var ii = 0; ii < OC_TUNE_ITRS; ii++) {
					const newTunedSine = tuneSineParams(bestSineParams, [MAIN_FREQUENCIES[s]*0.5, MAIN_FREQUENCIES[s]*2], OC_TUNE_TRIES/OC_TUNE_ITRS, chunkedXData, chunkedYData)
					if (newTunedSine.bestRss < bestRss) {
						bestRss = newTunedSine.bestRss
						bestSineParams = newTunedSine.params
					} else break // if there are no changes now, there will be no changes in the rest of the itrs
				}

				const thisFreq = bestSineParams[1] / 2 / Math.PI

				var thisMaximum = -bestSineParams[2] / bestSineParams[1] + Math.PI / bestSineParams[1] / 2
				if (thisMaximum < 0) thisMaximum += 2 * Math.PI / bestSineParams[1]
				maxima.push(thisMaximum + chunkXStart)
				// set back i so that it is a fourth of the period ahead of the found maximum
				while (dataToUse[0][i-1] > thisMaximum + chunkXStart + 1 / MAIN_FREQUENCIES[s] / 4) i--
				nextPeriodCutoff = thisMaximum + chunkXStart + 5 / MAIN_FREQUENCIES[s] / 4
				chunkedXData = []
				chunkedYData = []
				if (Date.now() >= nextUpdateTime) {
					console.log("Percent done", Math.round(i / dataToUse[0].length*100))
					nextUpdateTime += 2000
				}
			}
		}

		// claculate the oc values and make the graphs
		const ocValues = []
		for (var i = maxima.length-1; i >= 0; i--) {
			const oc = maxima[i] - maxima[i-1] - 1/MAIN_FREQUENCIES[s]
			if (Math.abs(oc) < 0.2) ocValues.unshift(oc*1440) // map OC to minutes
				else {
					maxima.splice(i, 1)
				}
		}

		console.log("The mean O-C is ", arrayMean(ocValues))
		console.log("The O-C stdev is ", arrayStdDev(ocValues, arrayMean(ocValues)))
		const lsrCoeffs = regression.linear(maxima.map((m, i) => [m, ocValues[i]]), { precision: 15 })
		console.log("there's a change of", lsrCoeffs.equation[0]*60, "seconds per day")
		console.log("the r2 is", lsrCoeffs.r2)
		const lsrPlotData = {
			x: [dataToUse[0][0], dataToUse[0][dataToUse[0].length-1]],
			y: [
				lsrCoeffs.predict(dataToUse[0][0])[1],
				lsrCoeffs.predict(dataToUse[0][dataToUse[0].length-1])[1]
			]
		}

		ocTraces.push({
			x: maxima.slice(1),
			y: ocValues,
			mode: "markers",
			marker: {
				color: "black",
				size: 4,
			},
			name: "OC",
			xaxis: "x"+String(s+1),
			yaxis: "y"+String(s+1)
		})

		bestFitTraces.push({
			x: lsrPlotData.x,
			y: lsrPlotData.y,
			mode: "lines",
			line: {
				color: "#e02b0f",
				width: 2.5
			},
			name: "best fit",
			xaxis: "x"+String(s+1),
			yaxis: "y"+String(s+1)
		})
	}

	const layout = {
		showlegend: false,
		margin: { t: 50 },
		xaxis: {
			title: 'Time (Days)',
			titlefont: {
				size: 20
			},
			tickfont: {
				size: 20
			}
		},
		yaxis: {
			title: 'O–C Period (Minutes)',
			titlefont: {
				size: 20
			},
			tickfont: {
				size: 20
			}
		},
		width: 704*2,
		height: 435,
		grid: {
			rows: 1,
			columns: 2,
			pattern: "independent"
		}
	}

	layout.xaxis2 = layout.xaxis
	layout.yaxis2 = layout.yaxis

	plot({
		data: [...ocTraces, ...bestFitTraces],
		layout,
		config: {staticPlot: true}
	})
}


function arrayMean(arr) {
	if (arr.length == 0) throw Error("Zero-length array")
	return arr.reduce((a, c) => a + c, 0) / arr.length
}

function arrayStdDev(arr, arrMean) {
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
		).slice(0, MAX_DATA_POINTS)

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

	return thisData
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
	return sineParams[0] * Math.sin(sineParams[1] * x + sineParams[2]) + sineParams[3]
}
function getRandomSineParams(rng, freqRange, dataRange) {
	const thisSineParams = []
	for (var i = 0; i < 4; i++) {
		const changeValue = rng.next()
		var newParam = 0
		if (i == 0) {
			newParam = 0.01 + (0.7 - 0.01) * changeValue
			newParam *= dataRange
		} else if (i == 1) { // variation by fixed interval for b (period)
			const thisFreq = freqRange[0] + (freqRange[1] - freqRange[0]) * changeValue
			newParam = 2 * Math.PI * thisFreq
		} else if (i == 2) { // phase angle is 0 to 2 pi
			newParam = 2 * Math.PI * changeValue
		}

		thisSineParams.push(newParam)
	}
	return thisSineParams
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

// tune the curve after random guesses
// tune count is per param
function tuneSineParams(initSineParams, freqRange, tuneCount, xData, yData) {
	const dataRange = Math.max(...yData) - Math.min(...yData)
	var bestRss = -1
	const bestSineParams = JSON.parse(JSON.stringify(initSineParams))
	for (var i = 0; i < 3; i++) {
		const thisSineParams = JSON.parse(JSON.stringify(bestSineParams))
		for (var ii = 0; ii < tuneCount; ii++) {
			const changeValue = ii / tuneCount
			var newParam = 0
			if (i == 0) {
				newParam = 0.01 + (0.7 - 0.01) * changeValue
				newParam *= dataRange
			} else if (i == 1) { // variation by fixed interval for b (period)
				const thisFreq = freqRange[0] + (freqRange[1] - freqRange[0]) * changeValue
				newParam = 2 * Math.PI * thisFreq
			} else if (i == 2) { // phase angle is 0 to 2 pi
				newParam = 2 * Math.PI * changeValue
			}

			thisSineParams[i] = newParam
			const thisRss = residualsSumOfSquares(thisSineParams, xData, yData)
			if (thisRss < bestRss || bestRss == -1) {
				bestRss = thisRss
				bestSineParams[i] = newParam
			}
		}
	}
	return { params: bestSineParams, bestRss: bestRss }
}

// find phase at certain JD time
function findPhaseAtJd(bParam, phaseAtZero, jd) {
	var phaseAtJd = bParam * jd + phaseAtZero
	while (phaseAtJd > 2 * Math.PI) {
		phaseAtJd -= 2 * Math.PI
	}
	return phaseAtJd
}

// find the line that goes through the data by binning
function findMiddleLine(binSize, inputX, inputY) {
	// find the bins
	const binX = [inputX[0]]
	const binY = [inputY[0]]
	var currentBinX = 0
	var currentBinY = 0
	var c = 0 // c for counter; how many elements are in bin already
	for (var i = 1; i < inputX.length; i++) {
		if (isNaN(inputX[i]) || isNaN(inputY[i])) {
			if (currentBinX != 0) {
				binX.push(currentBinX * binSize / c)
				binY.push(currentBinY * binSize / c)
				currentBinX = 0
				currentBinY = 0
			}
			binX.push(inputX[i-1])
			binY.push(NaN)
			c = 0
			continue
		}

		currentBinX += inputX[i] / binSize
		currentBinY += inputY[i] / binSize
		c ++
		if (c == binSize) {
			binX.push(currentBinX)
			binY.push(currentBinY)
			currentBinX = 0
			currentBinY = 0
			c = 0
		} else if (i == inputX.length-1) {
			binX.push(currentBinX * binSize / c)
			binY.push(currentBinY * binSize / c)
			currentBinX = 0
			currentBinY = 0
		}
	}

	binX.push(inputX[inputX.length - 1])
	binY.push(inputY[inputY.length - 1])

	return [binX, binY]
}

// find the noise
function findNoise(binX, binY, inputX, inputY) {
	const noiseY = []

	// calculate the noise
	var currentBinI = 0
	for (var i = 0; i < inputX.length; i++) {
		if (isNaN(inputX[i]) || isNaN(inputY[i])) {
			noiseY.push(NaN)
			throw Error("nan input")
		}

		if (inputX[i] > binX[currentBinI+1] && !realIsNaN(binY[currentBinI+2])) currentBinI ++
		if (inputX[i] > binX[currentBinI+2] && realIsNaN(binY[currentBinI+2])) currentBinI += 3
		const cb = currentBinI

		const lineSegSlope = (binY[cb+1] - binY[cb]) / (binX[cb+1] - binX[cb])
		const signalY = lineSegSlope * (inputX[i] - binX[cb]) + binY[cb]

		if (isNaN(signalY)) {
			console.log("nan signal y", i, cb, binX.length)
			console.log(binX[cb])
			continue
		}

		noiseY.push(inputY[i] - signalY)
	}

	return noiseY
}

function findNoisePower(noiseArr) {
	return noiseArr.reduce((a, c) => {
		if (isNaN(c)) return a
		return a + ((c**2) / noiseArr.length)
	}, 0)
}

// only checks if it is NaN
function realIsNaN(inp) {
	return isNaN(inp) && inp !== undefined && inp != null
}

function md5Hash(input) {
	return crypto.createHash("md5").update(input).digest("hex")
}
