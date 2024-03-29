import plotly as py
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def prepare_colors_and_lines(color_list):
    """
    Prepare a list of colors and a list of different line types for the plot.
    """

    colors = []
    lines = []

    line_types = ["solid", "dot", "dash", "longdash", "dashdot", "longdashdot"]

    for line in line_types:
        colors.extend(color_list)
        lines.extend([line] * len(color_list))

    return colors, lines

def plot_metagene_profiles(plot_df_start, plot_df_stop, read_length_list, title, max_y=None, color_list=None):
    """
    Create a scatter plot using plotly containing the sum of all metagene profiling lists split by alignment file for each chromosome.
    """

    x_range_start = plot_df_start["coordinates"]
    x_range_stop = plot_df_stop["coordinates"]

    labels = [label for label in plot_df_start.columns[1:].tolist() if int(label) in read_length_list]
    if max_y is None:
        max_y_start = plot_df_start[labels].max().max()
        max_y_stop = plot_df_stop[labels].max().max()
        max_y = max(max_y_start, max_y_stop)
        max_y += max_y * 0.1

    group_labels = [f"{count}" for count in range(len(labels))]

    fig = make_subplots(rows=1, cols=2, subplot_titles=("Start codon profile", "Stop codon profile"), horizontal_spacing=0.05)
    if color_list == []:
        color_list = py.colors.DEFAULT_PLOTLY_COLORS

    colors, lines = prepare_colors_and_lines(color_list)
    for i in range(len(labels)):
        fig.add_trace(go.Scatter(
            x=x_range_start,
            y=plot_df_start[labels[i]],
            mode='lines',
            name=labels[i],
            legendgroup=group_labels[i],
            line=dict(color=colors[i], dash=lines[i])
        ),
        row=1, col=1
    )
    fig.add_shape(
        go.layout.Shape(
            type="line",
            x0=0, y0=0, x1=0, y1=max_y,
            line = {"dash": "dash", "color": "grey", "width": 2}
        ),
        row=1, col=1
    )


    for i in range(len(labels)):
        fig.add_trace(go.Scatter(
            x=x_range_stop,
            y=plot_df_stop[labels[i]],
            mode='lines',
            name=labels[i],
            legendgroup=group_labels[i],
            showlegend=False,
            line=dict(color=colors[i], dash=lines[i])
        ),
        row=1, col=2
    )
    fig.add_shape(
        go.layout.Shape(
            type="line",
            x0=0, y0=0, x1=0, y1=max_y,
            line = {"dash": "dash", "color": "grey", "width": 2}
        ),
        row=1, col=2
    )


    fig.update_xaxes(
        title_text="Distance to start codon (nt)",
        title_font={"size": 20},
        row=1, col=1
    )

    fig.update_xaxes(
        title_text="Distance to stop codon (nt)",
        title_font={"size": 20},
        row=1, col=2
    )

    fig.update_yaxes(
        title_text="Read Coverage",
        title_font={"size": 20},
        range=[0, max_y],
        row=1, col=1
    )

    fig.update_yaxes(
        #title_text="Read Coverage",
        #title_font={"size": 20},
        range=[0, max_y],
        #side="right",
        row=1, col=2
    )

    fig.update_annotations(font_size=24)
    fig_layout = go.Layout(
        title={"text" : f"{title}", "x": 0.5, "xanchor": "center", "font": {"size": 28}}
    )
    fig.update_layout(fig_layout)

    return fig, max_y