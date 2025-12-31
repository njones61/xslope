# XSLOPE Input Template

The first step in using xslope is to load the input file. The input file contains all the necessary information about the slope, including the geometry, material properties, and loading conditions. The input file is an Excel file with a specific structure. A template for the Excel file can be downloaded here:

[input_template.xlsx](../../inputs/slope/input_template_lface5.xlsx)

The input file is designed to be easy to use and can be modified using any spreadsheet software. The input file is divided into several tabs or sheets, each of which corresponds to a specific aspect of the slope stability analysis. The strucutre of each sheet is designed to be intuitive and easy to understand. The sheets are as follows:

<div class="wrapped-table">
  <table>
    <thead>
      <tr>
        <th>Sheet Name</th>
        <th>Description</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>main</td>
        <td>A brief set of instructions and global variables including the unit weight of water, and tension crack parameters.</td>
      </tr>
      <tr>
        <td>plot</td>
        <td>A plot of the slope geometry based on the inputs in profile, piezo, and other sheets. This is intended to provide a quick visual check of the inputs.</td>
      </tr>
      <tr>
        <td>profile</td>
        <td>A set of tables for inputing the XY coordinates of up to 10 profile lines.</td>
      </tr>
      <tr>
        <td>mat</td>
        <td>A set of tables for inputting the material properties of up to 10 materials. This includes unit weight and shear strength properties.</td>
      </tr>
      <tr>
        <td>piezo</td>
        <td>A table for inputting the XY coordinates of piezometric line used to calculate pore pressures.</td>
      </tr>
      <tr>
        <td>circles</td>
        <td>A table for inputting the geometry of up to ten candidate circular failure surfaces.</td>
      </tr>
      <tr>
        <td>non-circ</td>
        <td>A table for inputting the XY coordinates of a non-circular slip surface.</td>
      </tr>
      <tr>
        <td>dloads</td>
        <td>A set of tables for inputting up to 8 distributed loads.</td>
      </tr>
      <tr>
        <td>reinforce</td>
        <td>A set of tables for inputting up to 12 soil reinforcement lines.</td>
      </tr>
    </tbody>
  </table>
</div>