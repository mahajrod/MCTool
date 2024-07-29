
# Main drawing objects (nested, from top to bottom)

1. Figure (Maybe remove)
2. Subplot
3. VerticalTrackGroup
4. HorizontalTrackGroup
5. Track

# Figure scheme

![image](figure_scheme.png)

1. Dimensions of the subplot:
   1. Horizontal dimensions
      1. S<sub>sl</sub>   - left space 
      2. S<sub>vlw</sub>  - width of the vertical track group label 
      3. S<sub>vls</sub>  - space between vertical track group label and connector between labels
      4. S<sub>cw</sub>   - connector width
      5. S<sub>cs</sub>   - space between connector and horizontal track group label
      6. S<sub>hlw</sub>  - width of horizontal track group label 
      7. S<sub>hls</sub>  - space between horizontal track group label and subplot's Y axis
      8. S<sub>ylw</sub>  - width of Y axis (including ticks and label)
      9. S<sub>ys</sub>   - space between Y axis and vertical track groups
      10. S<sub>vw</sub>  - width of vertical track groups
      11. S<sub>lsl</sub> - space between vertical track groups and legend
      12. S<sub>lw</sub>  - width of the legend
      13. S<sub>sr</sub>  - right space
   2. Vertical dimensions
      1. S<sub>vh1</sub>  - height of the first vertical track group
      2. S<sub>vs</sub>   - space between vertical track groups
      3. S<sub>vh2</sub>  - height of the second vertical track group
      
          ..................................................................
   
      4. S<sub>vhn</sub>  - height of the last vertical track group
      5. S<sub>ybs</sub>  - distance between left bottom corner of the last vertical track group and start of the subplot's Y axis
      6. S<sub>xst</sub>  - space between the last vertical track group and subplot's X axis 
      7. S<sub>xh</sub>   - height of the subplot's X axis (including ticks and label)
      8. S<sub>sb</sub>   - bottom space
2. Dimensions of the vertical track group
   1. Horizontal dimensions
      1. V<sub>sl</sub>   - left space 
      2. V<sub>yw</sub>   - width of the vertical track group's Y axis (including ticks and label)
      3. V<sub>ys</sub>   - space between vertical track group's Y axis and horizontal track groups
      4. V<sub>hw</sub>   - width of the horizontal track groups 
      5. V<sub>sr</sub>   - right space
      6. V<sub>w</sub>    - whole width of the vertical track groups
   2. Vertical dimensions
      1. V<sub>st</sub>   - top space
      2. V<sub>hh1</sub>  - height of the first horizontal track group
      3. V<sub>his</sub>  - space between horizontal track groups
      4. V<sub>hh2</sub>  - height of the second horizontal track group
    
         .....................................................................

      5. V<sub>Vhhn</sub>  - height of the last horizontal track group
      6. V<sub>Vxst</sub>  - space between the last horizontal track group and vertical track group's X axis 
      7. V<sub>Vxh</sub>   - height of the X axis (including ticks and label)
      8. V<sub>Vsb</sub>   - bottom space
3. Dimensions of the horizontal track group
   1. Horizontal dimensions
      1. H<sub>sl</sub>   - left space
      2. H<sub>yw</sub>   - width of the horizontal track group's Y axis (including ticks and label)
      3. H<sub>ys</sub>   - space between horizontal track group's Y axis and first track
      4. H<sub>tw1</sub>  - width of the first track 
      5. H<sub>sti</sub>  - space between tracks
      6. H<sub>tw2</sub>  - width of the second track
    
         .....................................................................

      7. H<sub>twn</sub>  - width of the last track
      8. H<sub>sb</sub>   - right space
   2. Vertical dimensions
      1. H<sub>st</sub>   - top space
      2. H<sub>ht</sub>   - height of the tracks
      3. H<sub>xst</sub>  - space between tracks and horizontal track group's X axis 
      4. H<sub>xh</sub>   - height of the horizontal track group's X axis (including ticks and label)
      5. H<sub>sb</sub>   - bottom space
4. Dimensions of the track
   1. Horizontal dimensions
      1. T<sub>sl</sub>   - left space
      2. T<sub>yw</sub>   - width of the track's Y axis (including ticks and label)
      3. T<sub>ys</sub>   - space between track's Y axis and track body
      4. T<sub>bw</sub>    - width of the track body
      5. T<sub>sr</sub>    - right space
   2. Vertical dimensions
      1. T<sub>st</sub>   - top space
      2. T<sub>lh</sub>   - height of the track label
      3. T<sub>lsb</sub>  - space between track label and track body
      4. T<sub>bh</sub>    - height of the  track body
      5. T<sub>xst</sub>  - space between track bodies and track's X axis 
      6. T<sub>xh</sub>   - height of the tracks' sX axis (including ticks and label)
      7. T<sub>sb</sub>   - bottom space

# Formulas:

*Track width* T<sub>w</sub> = T<sub>sl</sub> + T<sub>yw</sub> + T<sub>ys</sub> + T<sub>bw</sub> +  T<sub>sb</sub>

| dimension | T<sub>sl</sub>       | T<sub>yw</sub>  | T<sub>ys</sub>   | T<sub>bw</sub> | T<sub>sb</sub> |
|-----------|----------------------|-----------------|------------------|----------------|----------------|
| source    | style                | auto from style | style            | data           | style          |
| values    | fraction<sup>*</sup> | fraction        | fraction         | absolute       | fraction       |

\* fraction of max (T<sub>bw</sub>)

*Track height* T<sub>h</sub> = T<sub>st</sub> + T<sub>lh</sub> + T<sub>lsb</sub> + T<sub>bh</sub> + T<sub>xst</sub>  + T<sub>xh</sub> + T<sub>sb</sub>

| dimension | T<sub>st</sub> | T<sub>lh</sub>   | T<sub>lsb</sub> | T<sub>bh</sub> | T<sub>xst</sub> | T<sub>xh</sub>  | T<sub>sb</sub> |
|-----------|----------------|------------------|-----------------|----------------|-----------------|-----------------|----------------|
| source    | style          | auto from style  | style           | data           | style           | auto from style | style          |
| values    | absolute       | absolute         | absolute        | absolute       | absolute        | absolute        | absolute       |


*Horizontal track group width* H<sub>w</sub> = H<sub>sl</sub> + H<sub>yw</sub> + H<sub>ys</sub> + sum(T<sub>wi</sub>) + (N-1)* H<sub>sti</sub> + H<sub>sr</sub>

| dimension | H<sub>sl</sub>       | H<sub>yw</sub>  | H<sub>ys</sub> | T<sub>wi</sub> | H<sub>sti</sub> | H<sub>sr</sub> |
|-----------|----------------------|-----------------|----------------|----------------|-----------------|----------------|
| source    | style                | auto from style | style          | data           | style           | style          |
| values    | fraction<sup>*</sup> | fraction        | fraction       | absolute       | fraction        | fraction       |

\* fraction max(sum(T<sub>bw</sub>) per horizontal track group) 

*Horizontal track group height* H<sub>h</sub> = H<sub>st</sub> + H<sub>ht</sub> + H<sub>xst</sub> + H<sub>xh</sub> + H<sub>sb</sub>

| dimension | H<sub>st</sub> | H<sub>ht</sub>     | H<sub>xst</sub> | H<sub>xh</sub>  | H<sub>sb</sub> |
|-----------|----------------|--------------------|-----------------|-----------------|----------------|
| source    | style          | max(T<sub>h</sub>) | style           | auto from style | style          |
| values    | absolute       | absolute           | absolute        | absolute        | absolute       |

