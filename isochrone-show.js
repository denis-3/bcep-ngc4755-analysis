// show some isochrones from the astromancer files

const fs = require("fs")
const { plot } = require("nodeplotlib");

const FILE_PATHS = ["./cluster-data/our-isochrone.csv", "./cluster-data/lit-isochrone.csv"]

// Gaia filter
const BWCRU_XY = [0.148307, 9.031491]
const CVCRU_XY = [0.338291, 9.949834]
const COLOR_INDEX_NAME = 'BP–RP'
const LUM_FILTER_NAME = "G"

// BV filters
/*const BWCRU_XY = [0.06, 9.12]
const CVCRU_XY = [0.2, 10.49]
const COLOR_INDEX_NAME = 'B–V Magnitude'
const LUM_FILTER_NAME = "B"*/

const starPoints = []
const isochronePoints = []

FILE_PATHS.forEach((fn, i) => {
	const fData = fs.readFileSync(fn, "utf8").trim()
	fData.split("\n").slice(1).forEach(line => {
		const [sx, sy, ix, iy] = line.split(",")
		if (starPoints[i] === undefined) {
			starPoints[i] = {x: [], y: []}
			isochronePoints[i] = {x: [], y: []}
		}
		starPoints[i].x.push(Number(sx))
		starPoints[i].y.push(Number(sy))
		if (ix != "") {
			isochronePoints[i].x.push(Number(ix))
			isochronePoints[i].y.push(Number(iy))
		}
	})
})

var plotCount = 0
const traces = []
starPoints.forEach((s, i) => {
	traces.push({
		x: s.x,
		y: s.y,
		mode: "markers",
		marker: {
			color: "limegreen",
			opacity: 0.6,
			line: {
				color: "grey",
				width: 1
			}
		},
		xaxis: "x"+String(i+1),
		yaxis: "y"+String(i+1),
	})

	if (LUM_FILTER_NAME != "Bhh") {
		traces.push({
			x: isochronePoints[i].x,
			y: isochronePoints[i].y,
			mode: "lines",
			line: {
				color: "black",
				width: 2.5,
				shape: "spline",
				dash: "dot"
			},
			xaxis: "x"+String(i+1),
					yaxis: "y"+String(i+1)
		})
	}

	// custom star traces

	traces.push({
		x: [BWCRU_XY[0]],
		y: [BWCRU_XY[1]],
		mode: "markers",
		marker: {
			color: "#ff2a00",
			size: 12,
			symbol: "square",
			line: {
				color: "grey",
				width: 1
			}
		},
		xaxis: "x"+String(i+1),
		yaxis: "y"+String(i+1)
	})

	traces.push({
		x: [CVCRU_XY[0]],
		y: [CVCRU_XY[1]],
		mode: "markers",
		marker: {
			color: "#ff2a00",
			size: 14,
			symbol: "pentagon",
			line: {
				color: "grey",
				width: 1
			}
		},
		xaxis: "x"+String(i+1),
		yaxis: "y"+String(i+1)
	})

	plotCount ++
})

const rowCount = Math.ceil(traces.length / 8)

const layout = {
	showlegend: false,
	margin: {
		t: 50
	},
	xaxis: {
		title: `${COLOR_INDEX_NAME} Magnitude`,
		titlefont: {
			size: 20
		},
		tickfont: {
			size: 20
		}
	},
	yaxis: {
		autorange: "reversed",
		titlefont: {
			size: 20
		},
		tickfont: {
			size: 20
		}
	},
	width: plotCount > 1 ? 704*2 : 704,
	height: 435 * rowCount,
	grid: {
		rows: rowCount,
		columns: plotCount > 1 ? 2 : 1,
		pattern: "independent"
	}
}

for (var i = 2; i <= plotCount; i++) {
	layout[`xaxis${i}`] = JSON.parse(JSON.stringify(layout.xaxis))
	layout[`yaxis${i}`] = JSON.parse(JSON.stringify(layout.yaxis))
}

layout.yaxis.title = `${LUM_FILTER_NAME} Magnitude`

plot({layout, data: traces})
