const { isMainThread, workerData, parentPort } = require("worker_threads")

if (isMainThread) throw Error("Not supposed to be run as a main thread!")

var busy = false
const xData = new Float64Array(workerData.x)
const yData = new Float64Array(workerData.y)
const { freqRange, lombTries } = workerData

parentPort.on("message", (m) => {
	if (busy) throw Error("Prewhiten worker is busy!")
	busy = true
	if (m.command == "exit") {
		parentPort.removeAllListeners()
		process.exitCode = 0
		return
	}
	if (m.command != "start periodogram") throw Error("Unsupported command: " + m.command)
	var bestFreq = 0
	var bestPower = 0
	var periodogramX = []
	var periodogramY = []
	for (var ii = m.startI; ii < m.endI; ii++) {
		const thisFreq = freqRange[0] + (freqRange[1] - freqRange[0]) * ii / lombTries

		const thisTau = lombScargleTau(thisFreq, xData)
		const thisCosPart = lombScargleCosPart(thisFreq, xData, yData, thisTau)
		const thisSinPart = lombScargleSinPart(thisFreq, xData, yData, thisTau)
		const thisLombScarglePower = (thisCosPart + thisSinPart) / 2

		if (m.iteration == 0) periodogramX.push(thisFreq)
		periodogramY.push(2 * Math.sqrt(thisLombScarglePower / xData.length))
	}

	parentPort.postMessage({
		message: "finished periodogram",
		x: periodogramX,
		y: periodogramY,
	})
	busy = false
})

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
