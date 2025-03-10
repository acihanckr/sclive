import warnings
from typing import Optional, Union
import plotly.graph_objects as go
from ._layout_funcs import set_2d_layout
import polars as pl

def violin_plt(adata, 
            x_var:str, 
            meta_id:str, 
            group_by:Optional[str]=None, 
            layer:Optional[str]=None,
            vln_type:Optional[str]='single', 
            pts:Union[str,bool]=False, 
            pt_size:Optional[str]=4,
            jitter:Optional[float]=0.05,
            legend_size:Optional[str]=None,
            ticks_font_size:Optional[str]=12,
            width:Optional[str]='auto',
            height:Optional[str]='auto', 
            axis_font_size:Optional[str]=12)->go.Figure:
    '''
    Draws violin for a continuous observation meta or a gene expression for
    a given annotated data object

    Parameters:
    -----------
    :param adata: single cell object to be plotted
    :param x_var: x axis variable to draw violin/box plot
    :param meta_id: y axis variable to draw violin/box plot
    :param group_by: grouping variable for grouped violin/box plot
    :param vln_type: violin plot type. Options: 'single', 'grouped', 'split'
    :param pts: either to draw data points. Options: 'all', 'outliers', False
    :param pt_size: point size if data points are drawn. If None no points will be drawn.
    :param ticks_font_size: size of tick labels on x and y axis 
    :param txt_size: font size of the axis labels. If None, axis labels will be omitted 
    :param width: width of the plot. Can be auto or any value Plotly graph objects accepts
    :param height: height of the plot. Can be auto or any value Plotly graph objects accepts.
    ...
    Returns:
    --------
    plotly.graph_objects.Figure with desired violin
    '''

    if x_var not in adata.obs.columns.tolist():
        raise(ValueError("Given x variable is not found in Annotated Data!"))
    else:
        plotting_data = (pl.DataFrame(adata.obs.loc[:,[x_var]].reset_index()))
    
    if meta_id in adata.obs.columns.tolist():
        plotting_data = plotting_data.join(pl.DataFrame(adata.obs.loc[:,[meta_id]].reset_index()), on="index", how="inner")
    elif meta_id in adata.var_names.tolist():
        plotting_data = (plotting_data.join(
            pl.from_pandas(adata[:,meta_id].to_df(layer=layer).reset_index()), on="index", how="inner"))
    else:
        raise(ValueError("Given meta data or gene expression is not found in Annotated Data!"))
    
    if group_by is not None:
        if group_by not in adata.obs.columns.tolist():
            raise(ValueError("Given group by variable is not found in Annotated Data!"))
        plotting_data = (plotting_data.join(pl.DataFrame(adata.obs.loc[:,[group_by]].reset_index()), on="index", how="inner")
                         .rename({group_by:"group_by"}).with_columns(pl.col("group_by").cast(pl.String).cast(pl.Categorical)))
    elif vln_type != "single":        
        warnings.warn("Group by variable is not provided. Violin plot type will be set to single!")
        vln_type = "single"
    
    plotting_data = (plotting_data.rename({x_var: "X", meta_id: "Y", "index":"barcode"}, strict=False)
                        .with_columns(pl.col("X").cast(pl.String).cast(pl.Categorical),
                                        pl.col("Y").cast(pl.Float64)))

    fig = go.Figure()
    if vln_type=="single":
        for i in plotting_data["X"].unique():
            fig.add_trace(go.Violin(x=plotting_data.filter(pl.col("X") == i)["X"],
                            y=plotting_data.filter(pl.col("X") == i)["Y"],
                            marker=dict(size=pt_size),
                            points=pts, name=str(i)))
    elif vln_type=="grouped":
        for i in plotting_data["group_by"].unique():
            fig.add_trace(go.Violin(x=plotting_data.filter(pl.col("group_by") == i)["X"],
                            y=plotting_data.filter(pl.col("group_by") == i)["Y"],
                            marker=dict(size=pt_size),
                            points=pts, name=str(i)))
        fig.update_layout(violinmode="group")
    elif vln_type=="split":
        if len(plotting_data["group_by"].unique()) != 2:
            return fig
        else:
            groups = plotting_data["group_by"].unique().to_list()
            xs = plotting_data["X"].unique().to_list()
            show_legend = True
            for i in range(plotting_data.unique("X").shape[0]):
                fig.add_trace(go.Violin(x=plotting_data.filter((pl.col("group_by") == groups[0])&(pl.col("X") == xs[i]))["X"].to_list(),
                            y=plotting_data.filter((pl.col("group_by") == groups[0])&(pl.col("X") == xs[i]))["Y"].to_list(),
                            legendgroup=groups[0], scalegroup=groups[0], name=groups[0],
                            side='negative',
                            line_color='lightseagreen',
                            points=pts,
                            pointpos=-0.5,
                            jitter=jitter,
                            showlegend=show_legend
                ))
                fig.add_trace(go.Violin(x=plotting_data.filter((pl.col("group_by") == groups[0])&(pl.col("X") == xs[i]))["X"].to_list(),
                            y=plotting_data.filter((pl.col("group_by") == groups[1])&(pl.col("X") == xs[i]))["Y"].to_list(),
                            legendgroup=groups[1], scalegroup=groups[1], name=groups[1],
                            side='positive',
                            line_color='mediumpurple',
                            points=pts,
                            jitter=jitter,
                            pointpos=0.5,
                            showlegend=show_legend
                ))
                show_legend = False if show_legend else False
        fig.update_layout(violingap=0, violingroupgap=0, violinmode='overlay')
    
    fig = set_2d_layout(fig, 
                        ticks_font_size=ticks_font_size,
                        axis_font_size=axis_font_size,
                        title_size=None,
                        title=None,
                        dimred_labels=[x_var, meta_id],
                        legend_size=legend_size,
                        width=width,
                        height=height)
    
    return fig