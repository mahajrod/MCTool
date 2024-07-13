# General track format

MACE uses an extended BED format for storing tracks. 
It is a tab-separated format with first three columns containing id of scaffold, start and end of the 
region (zero-based coordinate system with half-open intervals, i,.e python-style notation.). The rest of columns contain track-specific values and modifiers.
A track file may contain multiple tracks of various types. 
First row of the file must start from '#' symbol and contain a header with column names.

**Header structure:**

```
#scaffold_id    start   end track1&type1 track1&par1 ... track1&parK ... trackN&typeN trackN&par1 ... trackN&parM
```
**Header fields**

| Column        | Description         | Values                                                   | Comment                                                                                                                                                                                                                       |
| ------------- |---------------------|----------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| scaffold_id | Id of the scaffold  | single str | NA                                                                                                                                                                                                                            |
| start | start of the region | single int | 0-based, start is included in the  region                                                                                                                                                                                     |
| end | end of the region   | single int | 0-based, end is not included in the region                                                                                                                                                                                    |
| track1&type1 | Track values for 'track1' | comma-separated list of values (str, float, int or bool) | Column name encodes type of the track. See **Allowed types of the tracks** for details. By default, *&* is used as separator between *track name* and *track type*                                                            |
| track1&par1 | Values for 'par1' of 'track1' | comma-separated list of values (str, float, int or bool) | Column name encodes parameter of the track. It should include *track name* and *parameter name*, separated by *&*. See **Recognizable track parameters** for details. If parameter is absent in the list, it will be ignored. |

**Track types:**

Some track types may include several internal separators (zero, one or more) "$" if it is necessary to set a specific subtype of the track or additional parameters.
Check specific tracks in **Allowed types of the tracks**. For example, "AAAAA&marker_m$ellipse" is a valid and meanins track type "marker_m" and marker type "ellipse" for track "AAAAA"  

**Attaching tracks:**

In some cases you may need to draw one track (query) over other one (target).
For this just add **@target_track_name** to the query track after the type.
If you wish to attach several query tracks to a single target track, then you will need  specify a query-target chain for each of queries in the order of drawing.  

*Example.* Two marker_m query tracks attached to the hist track. Query track with rectangle markers will be drawn last:

AAAAA$hist - target

BBBBB&marker_m$ellipse@AAAAA - query attached to the target AAAAA

CCCCC&marker_m$rectangle@BBBBB - query2 attached to the same target AAAAA, but will be drawn after query 1

**Allowed types of the tracks:**
1. **hist**          - use a full height of a track to plot a stacked histogram. Sum of components must not exceed 1.00 in any window. Use a script "prepare_hist_for_tracks.py" to prepare your values for the track.
2. **window**        - use a full height of the window to plot a rectangle.
3. **plot_m**        - draw a curve over whole track with points set at the middles of the regions
4. **plot_e**        - draw a curve over whole track with points set at the starts and ends of the regions
5. **marker_m**      - draw a marker at the middle of the region. Allows track subtypes. See **Allowed types of the markers**.
6. **marker_e**      - draw marker at each border of the region. Allows track subtypes. See **Allowed types of the markers**.

**Allowed types of the markers**:

1. **rect**      - draw a rectangle marker with start and stop on the edges of region, but with height smaller than track.
2. **ellipse**   - draw an ellipse with center at the middle between start and stop, height and x/y ratio are controlled by the feature style. Can be used to draw a circle if X and Y axis are in a different scale 
3. **circle**    - draw a circle with center at the middle between start and stop. It will look as an ellipse on the plot, if  X and Y axis are in a different scale. Use **ellipse** type to avoid this issue.

**Recognizable track parameters** implemented so far:
1. **colors** - comma-separated list of matplotlib color names or color hexs or keywords ("auto" or "default") 
2. **bg_color** - matplotlib color name or color hex or keyword ("default")
3. **height_offset** - height offset for the center of the marker. Fraction or keyword ("default")

**Keywords**:
1. **default** - use default from the script (vary from script to script)
2. **auto** - calculate according the logic implemented in the script (vary from script to script) if it exists. Otherwise, use default

**Example**:
```

```