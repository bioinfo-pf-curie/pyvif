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
