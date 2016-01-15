<table id="hor-minimalist-a">
% for key, val in info:
    <tr>
      <td><b>${key}</b></td>
      % if isinstance(val, list):
      <td>
        <table>
          % for v in val:
          <tr>
            <td>${v[0]}</td>
          </tr>
          % endfor
        </table>
      </td>
      % else:
      <td>${val}</td>
      % endif
    </tr>
% endfor
</table>
<pre>
% for seqcol in sequence:
    ${seqcol}
% endfor
</pre>
