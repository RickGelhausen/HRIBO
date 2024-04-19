#!/usr/bin/env python
import argparse
import pandas as pd

import plotly as py
import plotly.graph_objects as go

from pathlib import Path


INTRO_HTML = \
"""
<!DOCTYPE html>
<html>
    <head>
        <style>
            h1 { text-align: center; font-size: 3.5em; }
            h2 { text-align: center; font-size: 2em; }
            h3 { text-align: center; font-size: 1.5em; }
            h4 { text-align: left; font-size: 1.5em; }
            body { background-color: #f1f1f1; }
            .headerline { height: 2px; background-color: black; }
            .seperator { height: 1px; width: 90%;}
            .page { margin: 50px; }
            .description { margin-left: 100px; margin-right: 100px; margin-top: 0px; margin-bottom: 100px; }
            .plotly-graph-div { margin: 0 auto; }
        </style>
        <title>Differential Expression Quality Control</title>
    </head>
    <body>
        <div class="page">
"""

OUTRO_HTML = \
"""
        </div>
    </body>
</html>
"""


def custom_sort(group_by_names):
    """
    sort the group by names list based on condition
    """

    group_by_names = list(group_by_names)
    group_by_names.sort(key=lambda x: (x[0].split("_")[1], x[0].split("_")[0]), reverse=True)

    return group_by_names

def plot_scatter_2D(table_df, percentage_variance):
    """
    Plot the 2D scatter plot.
    """

    group_by_names = table_df.groupby("group")

    labels = []
    pca_data = []
    group_by_names = custom_sort(group_by_names)

    for group in group_by_names:
        labels.append(group[0])
        pca_data.append([group[1]["PC1"], group[1]["PC2"]])

    color_list = py.colors.DEFAULT_PLOTLY_COLORS

    group_labels = [f"{count}" for count in range(len(labels))]

    min_x = min([min(x[0]) for x in pca_data])
    max_x = max([max(x[0]) for x in pca_data])
    min_y = min([min(x[1]) for x in pca_data])
    max_y = max([max(x[1]) for x in pca_data])

    l_min_x = min_x + (min_x * 0.1) if (min_x * 0.1) < 0 else min_x - (min_x * 0.1)
    l_max_x = max_x - (max_x * 0.1) if (max_x * 0.1) < 0 else max_x + (max_x * 0.1)
    l_min_y = min_y + (min_y * 0.1) if (min_y * 0.1) < 0 else min_y - (min_y * 0.1)
    l_max_y = max_y - (max_y * 0.1) if (max_y * 0.1) < 0 else max_y + (max_y * 0.1)

    fig = go.Figure()
    for i in range(len(labels)):
        fig.add_trace(go.Scatter(
            x=pca_data[i][0],
            y=pca_data[i][1],
            mode="markers",
            marker=dict(color=color_list[i], size=12), #line=dict(color="rgb(0,0,0)", width=2)),
            name=labels[i],
            legendgroup=group_labels[i]
        ))

    fig.update_layout(
        xaxis = dict(#backgroundcolor="rgb(255,255,255)",
                    #gridcolor="rgb(0,0,0)",
                    title=f"PC1 ({percentage_variance[0] * 100:.2f}% variance)",
                    tickfont=dict(size=14),
                    range=[l_min_x, l_max_x],
                ),
        yaxis = dict(#backgroundcolor="rgb(255,255,255)",
                    #gridcolor="rgb(0,0,0)",
                    title=f"PC2 ({percentage_variance[1] * 100:.2f}% variance)",
                    tickfont=dict(size=14),
                    range=[l_min_y, l_max_y],
                ),
        font_family="Arial",
        font_size=16,
        legend=dict(
            itemsizing="constant",
            font=dict(size=18),
        ),
        margin=dict(t=50, b=20)
    )

    return fig

def plot_scatter_3D(table_df, percentage_variance):
    """
    Plot the 3D scatter plot.
    """
    group_by_names = table_df.groupby("group")

    labels = []
    pca_data = []
    group_by_names = custom_sort(group_by_names)

    for group in group_by_names:
        labels.append(group[0])
        pca_data.append([group[1]["PC1"], group[1]["PC2"], group[1]["PC3"]])

    color_list = py.colors.DEFAULT_PLOTLY_COLORS
    #color_list = ["rgb(0,0,0)", "rgb(183,183,183)", "rgb(250,128,114)", "rgb(238,238,0)", "rgb(65,105,225)", "rgb(189,252,201)"]

    group_labels = [f"{count}" for count in range(len(labels))]

    min_x = min([min(x[0]) for x in pca_data])
    max_x = max([max(x[0]) for x in pca_data])
    min_y = min([min(x[1]) for x in pca_data])
    max_y = max([max(x[1]) for x in pca_data])
    min_z = min([min(x[2]) for x in pca_data])
    max_z = max([max(x[2]) for x in pca_data])

    l_min_x = min_x + (min_x * 0.1) if (min_x * 0.1) < 0 else min_x - (min_x * 0.1)
    l_max_x = max_x - (max_x * 0.1) if (max_x * 0.1) < 0 else max_x + (max_x * 0.1)
    l_min_y = min_y + (min_y * 0.1) if (min_y * 0.1) < 0 else min_y - (min_y * 0.1)
    l_max_y = max_y - (max_y * 0.1) if (max_y * 0.1) < 0 else max_y + (max_y * 0.1)
    l_min_z = min_z + (min_z * 0.1) if (min_z * 0.1) < 0 else min_z - (min_z * 0.1)
    l_max_z = max_z - (max_z * 0.1) if (max_z * 0.1) < 0 else max_z + (max_z * 0.1)

    fig = go.Figure()
    for i in range(len(labels)):
        fig.add_trace(go.Scatter3d(
            x=pca_data[i][0],
            y=pca_data[i][1],
            z=pca_data[i][2],
            mode="markers",
            marker=dict(color=color_list[i]), #line=dict(color="rgb(0,0,0)", width=2)),
            name=labels[i],
            legendgroup=group_labels[i]
        ))

        for j in range(len(pca_data[i][0])):
            fig.add_trace(go.Scatter3d(
                x=[pca_data[i][0].tolist()[j], pca_data[i][0].tolist()[j]],
                y=[pca_data[i][1].tolist()[j], pca_data[i][1].tolist()[j]],
                z=[min_z, pca_data[i][2].tolist()[j]],
                mode="lines",
                showlegend=False,
                legendgroup=group_labels[i],
                line=dict(color="rgb(0,0,0)", width=2)
            ))

    fig.update_layout(
        scene = dict(
            xaxis = dict(#backgroundcolor="rgb(255,255,255)",
                        #gridcolor="rgb(0,0,0)",
                        title=f"PC1 ({percentage_variance[0] * 100:.2f} %)",
                        tickfont=dict(size=14),
                        range=[l_min_x, l_max_x],
                    ),
            yaxis = dict(#backgroundcolor="rgb(255,255,255)",
                        #gridcolor="rgb(0,0,0)",
                        title=f"PC2 ({percentage_variance[1] * 100:.2f} %)",
                        tickfont=dict(size=14),
                        range=[l_min_y, l_max_y],
                    ),
            zaxis = dict(#backgroundcolor="rgb(255,255,255)",
                        #gridcolor="rgb(0,0,0)",
                        title=f"PC3 ({percentage_variance[2] * 100:.2f} %)",
                        tickfont=dict(size=14),
                        range=[l_min_z, l_max_z],
                    ),
            camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
        ),
        font_family="Arial",
        font_size=16,
        legend=dict(
            itemsizing="constant",
            font=dict(size=18),
        )
    )

    return fig

def plot_correlation(cor_df):
    """
    Create a correlation heat map
    """

    fig = go.Figure()

    fig.add_trace(go.Heatmap(
        z=cor_df.values,
        x=[x.replace("_", "-") for x in list(cor_df.columns)],
        y=[y.replace("_", "-") for y in list(cor_df.index)],
        colorscale="Viridis",
        reversescale=True,
        zmax=1,
        zmin=0.5,
        text=cor_df.values,
        texttemplate="%{text:.2f}",
    ))

    fig.update_layout(
        xaxis=dict(side="top"),
        yaxis=dict(side="left", autorange="reversed"),
        font_family="Arial",
        font_size=16,
        margin=dict(t=50, b=20),
    )
    return fig

def create_html_file(fig_cor, fig_pca, output_path, file_suffix):
    """
    Write the html file.
    """
    html_string = INTRO_HTML + "\n"
    html_string += f"<h1>Differential expression Quality Control</h1>\n"
    html_string += f"<h2>Hierarchical Clustering Heatmap</h2>\n"
    html_string += fig_cor.to_html(config={"toImageButtonOptions": {"format" : "svg"}}, full_html=False, default_width="100%", default_height="800px" )
    html_string += "<div class=description>\n"
    html_string += f"<p>The heatmap displays the correlation of gene expression for all pairwise combinations of input samples. It indicates which samples are more similar to each other based on the normalized gene expression values.</p>\n"
    html_string += f"<p>Typically all samples have high correlations with each other (values >0.70). Samples that have a lower value may indicate an outlier in your data or sample contamination.</p>\n"
    html_string += "</div>\n"
    html_string += f"<h2>Principal Component Analysis (PCA)</h2>\n"
    html_string += fig_pca.to_html(config={"toImageButtonOptions": {"format" : "svg"}}, full_html=False, include_plotlyjs=False, default_width="100%", default_height="800px")
    html_string += "<div class=description>\n"
    html_string += f"<p>PCA is a statistical technique used to reduce the dimensionality of large data sets. It is does this by splitting the data into its principal components based on the variance in the data. <br>The first principal component (PC1) represents the maximum variance within the samples (PC2 the second highest etc...). PCA is a powerful tool to detect patterns and outliers among all the samples. </p>\n"
    html_string += "</div>\n"
    with open(output_path / f"diffex_QC{file_suffix}.html", "w") as f:
        f.write(html_string)

def create_html_file_3D(fig, output_path):
    """
    Write the html file.
    """

    full_html = fig.to_html(config={"toImageButtonOptions": {"format" : "svg"}} )
    with open(output_path / "PCA_3D.html", "w") as f:
        f.write(full_html)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Plot PCA as a 2D and 3D scatter plot.')

    parser.add_argument("-r", "--rld", action="store", dest="input_table", type=Path, required=True, help= "Table containing data for the 2D and 3D scatter plot.")
    parser.add_argument("-p", "--percentage_variance", action="store", dest="percentage_variance", type=Path, required=True, help= "Path input file containing the percentage variance per principle component.")
    parser.add_argument("-c", "--correlation_map", action="store", dest="correlation_map", type=Path, required=True, help= "Path input file containing the correlation map.")
    parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", type=Path, required=True, help= "Output directory for the scatter plots file.")
    parser.add_argument("-s", "--file_suffix", action="store", dest="file_suffix", type=str, default="", help= "Suffix to add to the output file name.")
    args = parser.parse_args()

    # read the table
    table_df = pd.read_csv(args.input_table, sep="\t")

    cor_df = pd.read_csv(args.correlation_map, sep="\t")
    fig_cor = plot_correlation(cor_df)

    with open(args.percentage_variance, "r") as f:
        percentage_variance = [float(line.rstrip()) for line in f]

    fig_2D = plot_scatter_2D(table_df, percentage_variance)
    fig_3D = plot_scatter_3D(table_df, percentage_variance)

    create_html_file(fig_cor, fig_2D, args.output_dir, f"{args.file_suffix}")
    create_html_file_3D(fig_3D, args.output_dir)

if __name__ == '__main__':
    main()