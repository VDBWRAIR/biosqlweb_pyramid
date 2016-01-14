<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html>
  <head>
    <link href="${request.static_url('biosql_pyramid:static/css/biosqlweb.css')}" rel="stylesheet">
    <link href="${request.static_url('biosql_pyramid:static/css/jquery/jquery-ui-1.7.custom.css')}" rel="stylesheet">

    <script src="${request.static_url('biosql_pyramid:static/javascript/jquery/jquery-1.3.2.min.js')}"></script>
    <script src="${request.static_url('biosql_pyramid:static/javascript/jquery/jquery-ui-1.7.custom.min.js')}"></script>
    <script src="${request.static_url('biosql_pyramid:static/javascript/jquery/ajaxupload.2.1.js')}"></script>
    <script src="${request.static_url('biosql_pyramid:static/javascript/summary.js')}"></script>
  <title>BioSQL Web</title>
</head>
<body>
  <h1>BioEntry Records</h1>
  <div id="bioentry_records">
    % for record in records:
      ${record | n}
    % endfor
  </div>
  <br/>
  <a href="#" id="button1" class="fg-button ui-state-default fg-button-icon-left ui-corner-all"><span class="ui-icon ui-icon-newwin"></span>GenBank Upload</a>
</body>
</html>
