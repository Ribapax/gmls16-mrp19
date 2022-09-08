const ChartJSImage = require("chart.js-image");


const plotData = (
    dataSetV1,
    dataSetV2,
    title,
    unitOfMeasurement
) => {
    return ChartJSImage().chart({
        "type": "line",
        "data": {
            "labels": ["32", "33", "64", "65", "128", "129", "256", "257", "512", "1000", "2000", "4000", "6000", "10000"],
            "datasets": [
                {
                    "label": "Antes da otimização",
                    "borderColor": "rgb(255,+99,+132)",
                    "backgroundColor": "rgba(255,+99,+132,+.5)",
                    "data": dataSetV1
                },
                {
                    "label": "Depois da otimização",
                    "borderColor": "rgb(54,+162,+235)",
                    "backgroundColor": "rgba(54,+162,+235,+.5)",
                    "data": dataSetV2
                },
            ]
        },
        "options": {
            "title": {
                "display": true,
                "text": title
            },
            "scales": {
                "xAxes": [
                    {
                        "display": true,
                        "scaleLabel": {
                            "display": true,
                            "labelString": "Tamanho da matriz"
                        }
                    }
                ],
                "yAxes": [
                    {
                        "display": true,
                        "type": "logarithmic",
                        "stacked": true,
                        "scaleLabel": {
                            "display": true,
                            "labelString": unitOfMeasurement
                        }
                    }
                ]
            }
        }
    }).backgroundColor('white')
        .width(500)
        .height(300);
}


async function main () {

    let v1DataSet = [57, 90, 11, 15, 37, 37, 27, 57, 90, 15, 11, 37, 37, 27]
    let v2DataSet = [71, 36, 94, 78, 98, 65, 61, 71, 36, 94, 78, 98, 65, 61]
    let title = "Sistema Linear - Grupo FLOPS_DP - Indicador DP MFLOPS/s"
    let unitOfMeasurement = "MFLOPS/s"

    const line_chart = plotData(
        v1DataSet,
        v2DataSet,
        title,
        unitOfMeasurement
    )

    await line_chart.toFile('chart.png');
}
main()