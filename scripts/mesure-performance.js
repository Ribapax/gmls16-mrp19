
const { promisify } = require('util')
const { parse } = require("csv-parse");
const fs = require('fs');

const exec = promisify(require('child_process').exec)

const csvMock = ``

const ITERATIONS_LIMIT = 10

const LIKWID_COMMAND = 'likwid-perfctr'
const FIRST_FLAGS = '-O -C 1 -g'
const SECOND_FLAGS = '-m'
const PROGRAM = 'invmat'

const TIME = 'TIME'
const L3 = 'L3'
const L2CACHE = 'L2CACHE'
const FLOPS_DP = 'FLOPS_DP'
const AVX_FLOPS_DP = 'AVX_FLOPS_DP'

const groups = [TIME, L3, L2CACHE, FLOPS_DP, AVX_FLOPS_DP]
const sizes = [32, 33, 64, 65, 128, 129, 256, 257, 512, 1000, 2000, 4000, 6000, 10000]

const getFileContents = async (filepath) => {
    const data = [];
    return new Promise(function(resolve, reject) {
        fs.createReadStream(filepath)
            .pipe(parse({delimiter: ',', relax_column_count: true, relax_quotes: true}))
            .on('error', error => reject(error))
            .on('data', row => data.push(row))
            .on('end', () => {
                resolve(data);
            });
    });
}

const parseKey = async (key, file) => {

    const data = await getFileContents(file)

    if (key === "TIME") {
        let lsTime, resTime = 0;
        for (let i = 0; i < data.length; i++) {
            if (data[i][0].search('Tempo iter') !== -1) {
                lsTime = +(data[i][0].substring(data[i][0].indexOf(':')+1))
                resTime = +(data[i+1][0].substring(data[i][0].indexOf(':')+1))
                break;
            }
        }
        return {
            'linearSystem': lsTime,
            'residue': resTime
        }
    }

    let linearSystemCalculationL3 = 0;
    let i;
    for (i = 0; i < data.length; i++) {
        if (data[i][0] === key) {
            linearSystemCalculationL3 = +(data[i][1])
            break;
        }
    }

    let residueL3 = 0;
    for (; i < data.length; i++) {
        if (data[i][0] === key) {
            residueL3 = +(data[i][1])
        }
    }

    return {
        'linearSystem': linearSystemCalculationL3,
        'residue': residueL3
    }
}

const parseTime = async (size) => {
    return parseKey('TIME', `output-${TIME}-${size}.csv`)
}

const parseL3 = async (size) => {
    return parseKey('L3 bandwidth [MBytes/s]', `output-${L3}-${size}.csv`)
}

const parseL2Cache = (size) => {
    return parseKey('L2 miss ratio', `output-${L2CACHE}-${size}.csv`)
}

const parseFlopsDP = (size) => {
    return parseKey('DP MFLOP/s', `output-${FLOPS_DP}-${size}.csv`)
}

const parseAVXFlopsDP = (size) => {
    return parseKey('AVX DP MFLOP/s', `output-${AVX_FLOPS_DP}-${size}.csv`)
}

const parsers = {
    L3: parseL3,
    L2CACHE: parseL2Cache,
    FLOPS_DP: parseFlopsDP,
    AVX_FLOPS_DP: parseAVXFlopsDP,
    TIME: parseTime,
}

const buildCommand = (group, size) => {
    let outputFileName = `output-${group}-${size}.csv`
    group = group === 'AVX_FLOPS_DP' || "TIME" ? 'FLOPS_DP' : group // Technical Resource
    return `${LIKWID_COMMAND} ${FIRST_FLAGS} ${group} ${SECOND_FLAGS} ./${PROGRAM} -r ${size} -i ${ITERATIONS_LIMIT} -s invmat-output > ${outputFileName}`
}

const execMock = async (command) => {
    console.log('Executing command: ' + command)
    return {
        stdout: csvMock
    }
}

const run = async (group, size, parser, mockExecution) => {
    const command = buildCommand(group, size)
    try {
        if (mockExecution) {
            await execMock(command)
            return parser(size)
        }
        console.log('Executing command: ' + command)
        await exec(command)
        return parser(size)
    } catch (e) {
        console.error(e)
        return parser(size)
    }
}

function Result(size, indicators) {
    this.matrixSize = size
    for (const indicator of indicators) {
        this[indicator.group] = indicator.value
    }
}

const main = async () => {

    const mockExecution = process.argv[2] === 'mock'
    const fixedSize = parseInt(process.argv[3])
    let actualSizes = sizes;
    if (fixedSize > 0) {
        actualSizes = [fixedSize]
    }

    const resultsPromises = actualSizes.map(async (size) => {

        const indicatorsPromises = groups.map(async group => ({
            size: size,
            group: group,
            values: await run(group, size, parsers[group], mockExecution)
        }))

        let indicators = []
        let i = 0;
        for (const promise of indicatorsPromises) {
            indicators[i] = await Promise.resolve(promise)
            i++
        }

        const linearSystemIndicators = indicators.map(indicator => {
            return {
                size: indicator.size,
                group: indicator.group,
                value: indicator.values.linearSystem
            }
        })

        const residueIndicators = indicators.map(indicator => {
            return {
                size: indicator.size,
                group: indicator.group,
                value: indicator.values.residue
            }
        })

        return {
            linearSystemResult: new Result(size, linearSystemIndicators),
            residueResult: new Result(size, residueIndicators),
        }
    })

    let results = []
    let i = 0;
    for (const promise of resultsPromises) {
        results[i] = await Promise.resolve(promise)
        console.log('\nLINEAR SYSTEM PARTIAL RESULT:\n')
        console.table(results[i].linearSystemResult)
        console.log('\nRESIDUE PARTIAL RESULT:\n')
        console.table(results[i].residueResult)
        i++
    }

    console.log('\n============== FINAL RESULTS ==================\n')

    const linearSystemResults = results.map(result => result.linearSystemResult)
    const residueResults = results.map(result => result.residueResult)

    console.log('\nLINEAR SYSTEM RESULTS:\n')

    console.table(linearSystemResults)

    console.log('\nRESIDUE RESULTS:\n')

    console.table(residueResults)

    // Output dataset for plot

    groups.forEach(group => {
        console.log("GROUP: " + group)
        const arr = linearSystemResults.map(result => result[group])
        console.log("Linear System: ")
        console.log(arr)

        const arr2 = residueResults.map(result => result[group])
        console.log("Residue: ")
        console.log(arr2)
        console.log("\n")
    })

}

main()