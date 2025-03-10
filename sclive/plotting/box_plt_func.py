import warnings
from typing import Optional, Union
import plotly.graph_objects as go
from ._layout_funcs import set_2d_layout
import polars as pl

def box_plt(adata, 
            x_var:str, 
            meta_id:str, 
            group_by:Optional[str]=None, 
            layer:Optional[str]=None,
            box_type:Optional[str]=None, 
            pts:Union[str,bool]=False, 
            pt_size:Optional[str]=4,
            jitter:Optional[float]=0.05,
            legend_size:Optional[str]=None,
            ticks_font_size:Optional[str]=12,
            width:Optional[str]='auto',
            height:Optional[str]='auto', 
            axis_font_size:Optional[str]=12)->go.Figure:
    '''
    Draws boxplot for a continuous observation meta or a gene expression for
    a given annotated data object

    Parameters:
    -----------
    :param adata: single cell object to be plotted
    :param x_var: x axis variable to draw violin/box plot
    :param meta_id: y axis variable to draw violin/box plot
    :param group_by: grouping variable for grouped violin/box plot
    :param box_type: box plot type. Options: 'single', 'grouped'
    :param pts: either to draw data points. Options: 'all', 'outliers', False
    :param pt_size: point size if data points are drawn. If None no points will be drawn.
    :param ticks_font_size: size of tick labels on x and y axis 
    :param txt_size: font size of the axis labels. If None, axis labels will be omitted 
    :param width: width of the plot. Can be auto or any value Plotly graph objects accepts
    :param height: height of the plot. Can be auto or any value Plotly graph objects accepts.
    ...
    Returns:
    --------
    plotly.graph_objects.Figure with desired boxplot
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
    elif box_type != "single":
        warnings.warn("Group by variable is not provided. Box plot type will be set to single!")
        box_type = "single"
    
    plotting_data = (plotting_data.rename({x_var: "X", meta_id: "Y", "index":"barcode"}, strict=False)
                        .with_columns(pl.col("X").cast(pl.String).cast(pl.Categorical),
                                        pl.col("Y").cast(pl.Float64)))

    fig = go.Figure()
    if box_type=="single":
        for i in plotting_data["X"].unique():
            fig.add_trace(go.Box(x=plotting_data.filter(pl.col("X") == i)["X"],
                            y=plotting_data.filter(pl.col("X") == i)["Y"],
                                name=str(i), marker=dict(size=pt_size),
                                boxpoints=pts))
    elif box_type=="grouped":
        for i in plotting_data["group_by"].unique():
            fig.add_trace(go.Box(x=plotting_data.filter(pl.col("group_by") == i)["X"],
                            y=plotting_data.filter(pl.col("group_by") == i)["Y"],
                                name=str(i), marker=dict(size=pt_size),
                                boxpoints=pts))
            fig.update_layout(boxmode="group")
    
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