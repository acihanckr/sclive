import warnings
from typing import Optional, Union
import plotly.graph_objects as go
from ._layout_funcs import set_2d_layout
import polars as pl
from sclive.dataio.get_metas_func import get_metas
from sclive.dataio.get_gene_exprs_func import get_gene_exprs

def box_plt(adata, 
            x_var:str, 
            meta_id:str, 
            group_by:Optional[str]=None, 
            layer:Optional[str]=None,
            use_raw:Optional[bool] = None,
            box_type:Optional[str]=None, 
            pts:Union[str,bool]=False, 
            pt_size:Optional[str]=4,
            legend_size:Optional[str]=None,
            ticks_font_size:Optional[str]=12,
            title:Optional[str]=None,
            title_size:Optional[int] = None,
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

    if meta_id in adata.obs.columns.tolist():
        plotting_data = get_metas(adata, [x_var], True).join(get_metas(adata, [meta_id], False), on="barcode")
    elif meta_id in adata.var_names.tolist():
        plotting_data = get_metas(adata, [x_var], True).join(get_gene_exprs(adata, [meta_id], layer=layer, use_raw=use_raw), on="barcode")
    else:
        raise(ValueError("Given meta data or gene expression is not found in Annotated Data!"))
    
    if group_by is not None:
        if group_by not in adata.obs.columns.tolist():
            raise(ValueError("Given group by variable is not found in Annotated Data!"))
        plotting_data = plotting_data.join(get_metas(adata, [group_by], True), on="barcode", how="inner")
    elif box_type != "single":
        warnings.warn("Group by variable is not provided. Box plot type will be set to single!")
        box_type = "single"
    
    fig = go.Figure()
    if box_type=="single":
        for i in plotting_data[x_var].unique():
            fig.add_trace(go.Box(x=plotting_data.filter(pl.col(x_var) == i)[x_var],
                            y=plotting_data.filter(pl.col(x_var) == i)[meta_id],
                                name=str(i), marker=dict(size=pt_size),
                                boxpoints=pts))
    elif box_type=="grouped":
        for i in plotting_data[group_by].unique():
            fig.add_trace(go.Box(x=plotting_data.filter(pl.col(group_by) == i)[x_var],
                            y=plotting_data.filter(pl.col(group_by) == i)[meta_id],
                                name=str(i), marker=dict(size=pt_size),
                                boxpoints=pts))
            fig.update_layout(boxmode="group")
    
    if title_size is not None and title is None:
        title = f'{meta_id} vs {x_var} grouped by {group_by} Box Plot' if group_by else f'{meta_id} vs {x_var} Box Plot'
    fig = set_2d_layout(fig, 
                        ticks_font_size=ticks_font_size,
                        axis_font_size=axis_font_size,
                        title_size=title_size,
                        title=title,
                        dimred_labels=[x_var, meta_id],
                        legend_size=legend_size,
                        width=width,
                        height=height)
    
    return fig