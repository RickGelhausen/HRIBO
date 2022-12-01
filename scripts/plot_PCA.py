#!/usr/bin/env python
import argparse
import pandas as pd

import plotly as py
import plotly.graph_objects as go

from pathlib import Path

def custom_sort(group_by_names):
    """
    sort the group by names list based on condition
    """

    group_by_names = list(group_by_names)
    group_by_names.sort(key=lambda x: (x[0].split("_")[1], x[0].split("_")[0]), reverse=True)

    return group_by_names

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

def create_html_file(fig, output_path):
    """
    Write the html file.
    """

    full_html = fig.to_html(config={"toImageButtonOptions": {"format" : "svg"}} )
    with open(output_path / "PCA_3D.html", "w") as f:
        f.write(full_html)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Plot PCA as a 3D scatter plot.')

    parser.add_argument("-r", "--rld", action="store", dest="input_table", type=Path, required=True, help= "Table containing data for the 3D scatter plot.")
    parser.add_argument("-p", "--percentage_variance", action="store", dest="percentage_variance", type=Path, required=True, help= "Path input file containing the percentage variance per principle component.")
    parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", type=Path, required=True, help= "Output directory for the scatter plots file.")
    args = parser.parse_args()

    # read the table
    table_df = pd.read_csv(args.input_table, sep="\t")

    with open(args.percentage_variance, "r") as f:
        percentage_variance = [float(line.rstrip()) for line in f]

    fig = plot_scatter_3D(table_df, percentage_variance)
    # fig.write_image(args.output_dir / "PCA_3D.pdf")
    # fig.write_image(args.output_dir / "PCA_3D.svg")
    create_html_file(fig, args.output_dir)


if __name__ == '__main__':
    main()
