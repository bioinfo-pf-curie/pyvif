// JS for pyvif report

$(function () {
  // enable the bootstrap tooltip hovers
  $('[data-toggle="tooltip"]').tooltip();

  // basic function to hide side navigation-bar.
  $('#sidenav-collapse').on('click', function () {
    $('#sidenav, #mainpage, #sidenav-logo').toggleClass('hide-nav');
    $('#sidenav-collapse span').toggleClass('fa-angle-left fa-angle-right');
  });
});

$(document).ready(function(){
  $('#bp_clustering').DataTable({
    'dom': '<"top"Bfi>rt',
    'order': [[ 6, 'desc' ]],
    'buttons': ['csv', 'copy'],
    'deferRender': true,
    'scrollY': 300,
    'scrollCollapse': true,
    'scroller': true
  });
});
