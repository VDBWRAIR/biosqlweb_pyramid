<script src="${request.static_url('biosql_pyramid:static/javascript/jquery/jquery-1.3.2.min.js')}"></script>
<div id='bioentry-${id}'>
  <div style='float: right'>
    <button id="download-fasta-${id}">Download Fasta</button>
  </div>
  <div>
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
  </div>
</div>
<script type='text/javascript'>
  // Provide download fasta interaction
  $( ":button" ).click(function() {
    accession = this.id.split('-')[2];
    window.location.href = 'download_file?accession=' + accession + '&format=fasta'
  });
</script>
