import pandas as pd
import plotly.graph_objects as go

# Example data (replace with your own)
# Suppose your DataFrame has columns: 'time', 'x', 'y'
# Suppose your DataFrame has columns: 'time', 'x', 'y'
df = pd.read_csv("x_2.csv")

# Create one scatter per time step (for slider)
frames = [
    go.Frame(
        data=[go.Scatter(x=[df.loc[i, 'x']], y=[df.loc[i, 'y']], mode='markers')],
        name=str(df.loc[i, 'time'])
    )
    for i in range(len(df))
]

# Initial plot
fig = go.Figure(
    data=[go.Scatter(x=[df.loc[0, 'x']], y=[df.loc[0, 'y']], mode='markers')],
    layout=go.Layout(
        xaxis=dict(range=[df['x'].min() - 1, df['x'].max() + 1]),
        yaxis=dict(range=[df['y'].min() - 1, df['y'].max() + 1]),
        updatemenus=[dict(
            type='buttons',
            showactive=False,
            buttons=[dict(label='Play',
                          method='animate',
                          args=[None, {"frame": {"duration": 500, "redraw": True},
                                       "fromcurrent": True}])]
        )]
    ),
    frames=frames
)

# Add slider
fig.update_layout(
    sliders=[{
        "steps": [{
            "args": [[str(t)], {"frame": {"duration": 0, "redraw": True}}],
            "label": str(t),
            "method": "animate"
        } for t in df['time']],
        "transition": {"duration": 0},
        "x": 0.1,
        "len": 0.9
    }]
)

fig.show()