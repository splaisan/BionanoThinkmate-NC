###########################################################################
#
# This file is auto-generated by the Perl DateTime Suite locale
# generator (0.05).  This code generator comes with the
# DateTime::Locale distribution in the tools/ directory, and is called
# generate-from-cldr.
#
# This file as generated from the CLDR XML locale data.  See the
# LICENSE.cldr file included in this distribution for license details.
#
# This file was generated from the source file si.xml
# The source file version number was 1.7, generated on
# 2009/05/05 23:06:40.
#
# Do not edit this file directly.
#
###########################################################################

package DateTime::Locale::si;

use strict;
use warnings;
use utf8;

use base 'DateTime::Locale::root';

sub cldr_version { return "1\.7\.1" }

{
    my $am_pm_abbreviated = [ "පෙ\.ව\.", "ප\.ව\." ];
    sub am_pm_abbreviated { return $am_pm_abbreviated }
}
{
    my $date_format_full = "EEEE\,\ y\ MMMM\ d";
    sub date_format_full { return $date_format_full }
}

{
    my $date_format_long = "y\ MMMM\ d";
    sub date_format_long { return $date_format_long }
}

{
    my $date_format_medium = "y\ MMM\ d";
    sub date_format_medium { return $date_format_medium }
}

{
    my $date_format_short = "yyyy\/MM\/dd";
    sub date_format_short { return $date_format_short }
}

{
    my $day_format_abbreviated = [ "සඳු", "අඟ", "බදා", "බ්‍රහ", "සිකු", "සෙන", "ඉරි" ];
    sub day_format_abbreviated { return $day_format_abbreviated }
}

sub day_format_narrow { $_[0]->day_stand_alone_narrow() }

{
    my $day_format_wide = [ "සඳුදා", "අඟහරුවාදා", "බදාදා", "බ්‍රහස්පතින්දා", "සිකුරාදා", "සෙනසුරාදා", "ඉරිදා" ];
    sub day_format_wide { return $day_format_wide }
}

sub day_stand_alone_abbreviated { $_[0]->day_format_abbreviated() }

{
    my $day_stand_alone_narrow = [ "ස", "අ", "බ", "බ්‍ර", "සි", "සෙ", "ඉ" ];
    sub day_stand_alone_narrow { return $day_stand_alone_narrow }
}

sub day_stand_alone_wide { $_[0]->day_format_wide() }

{
    my $era_abbreviated = [ "ක්‍රි\.පූ\.", "ක්‍රි\.ව\." ];
    sub era_abbreviated { return $era_abbreviated }
}

sub era_narrow { $_[0]->era_abbreviated() }

{
    my $era_wide = [ "ක්‍රිස්තු\ පූර්‍ව", "ක්‍රිස්තු\ වර්‍ෂ" ];
    sub era_wide { return $era_wide }
}
{
    my $first_day_of_week = "1";
    sub first_day_of_week { return $first_day_of_week }
}

{
    my $month_format_abbreviated = [ "ජන", "පෙබ", "මාර්ත", "අප්‍රේල", "මැය", "ජූන", "ජූල", "අගෝ", "සැප", "ඔක", "නොවැ", "දෙසැ" ];
    sub month_format_abbreviated { return $month_format_abbreviated }
}

sub month_format_narrow { $_[0]->month_stand_alone_narrow() }

{
    my $month_format_wide = [ "ජනවාර", "පෙබරවාර", "මාර්ත", "අප්‍රේල්", "මැයි", "ජූන", "ජූලි", "අගෝස්තු", "සැප්තැම්බර්", "ඔක්තෝබර්", "නොවැම්බර්", "දෙසැම්බර්" ];
    sub month_format_wide { return $month_format_wide }
}

sub month_stand_alone_abbreviated { $_[0]->month_format_abbreviated() }

{
    my $month_stand_alone_narrow = [ "ජ", "පෙ", "මා", "අ", "මැ", "ජූ", "ජූ", "අ", "සැ", "ඔ", "නො", "දෙ" ];
    sub month_stand_alone_narrow { return $month_stand_alone_narrow }
}

sub month_stand_alone_wide { $_[0]->month_format_wide() }

{
    my $quarter_format_abbreviated = [ "කාර්\:1", "කාර්\:2", "කාර්\:3", "කාර්\:4" ];
    sub quarter_format_abbreviated { return $quarter_format_abbreviated }
}
{
    my $quarter_format_wide = [ "1\ වන\ කාර්තුව", "2\ වන\ කාර්තුව", "3\ වන\ කාර්තුව", "4\ වන\ කාර්තුව" ];
    sub quarter_format_wide { return $quarter_format_wide }
}

sub quarter_stand_alone_abbreviated { $_[0]->quarter_format_abbreviated() }


sub quarter_stand_alone_wide { $_[0]->quarter_format_wide() }

{
    my $time_format_full = "h\:mm\:ss\ a\ zzzz";
    sub time_format_full { return $time_format_full }
}

{
    my $time_format_long = "h\:mm\:ss\ a\ z";
    sub time_format_long { return $time_format_long }
}

{
    my $time_format_medium = "h\:mm\:ss\ a";
    sub time_format_medium { return $time_format_medium }
}

{
    my $time_format_short = "h\:mm\ a";
    sub time_format_short { return $time_format_short }
}

1;

__END__


=pod

=encoding utf8

=head1 NAME

DateTime::Locale::si

=head1 SYNOPSIS

  use DateTime;

  my $dt = DateTime->now( locale => 'si' );
  print $dt->month_name();

=head1 DESCRIPTION

This is the DateTime locale package for Sinhala.

=head1 DATA

This locale inherits from the L<DateTime::Locale::root> locale.

It contains the following data.

=head2 Days

=head3 Wide (format)

  සඳුදා
  අඟහරුවාදා
  බදාදා
  බ්‍රහස්පතින්දා
  සිකුරාදා
  සෙනසුරාදා
  ඉරිදා

=head3 Abbreviated (format)

  සඳු
  අඟ
  බදා
  බ්‍රහ
  සිකු
  සෙන
  ඉරි

=head3 Narrow (format)

  ස
  අ
  බ
  බ්‍ර
  සි
  සෙ
  ඉ

=head3 Wide (stand-alone)

  සඳුදා
  අඟහරුවාදා
  බදාදා
  බ්‍රහස්පතින්දා
  සිකුරාදා
  සෙනසුරාදා
  ඉරිදා

=head3 Abbreviated (stand-alone)

  සඳු
  අඟ
  බදා
  බ්‍රහ
  සිකු
  සෙන
  ඉරි

=head3 Narrow (stand-alone)

  ස
  අ
  බ
  බ්‍ර
  සි
  සෙ
  ඉ

=head2 Months

=head3 Wide (format)

  ජනවාර
  පෙබරවාර
  මාර්ත
  අප්‍රේල්
  මැයි
  ජූන
  ජූලි
  අගෝස්තු
  සැප්තැම්බර්
  ඔක්තෝබර්
  නොවැම්බර්
  දෙසැම්බර්

=head3 Abbreviated (format)

  ජන
  පෙබ
  මාර්ත
  අප්‍රේල
  මැය
  ජූන
  ජූල
  අගෝ
  සැප
  ඔක
  නොවැ
  දෙසැ

=head3 Narrow (format)

  ජ
  පෙ
  මා
  අ
  මැ
  ජූ
  ජූ
  අ
  සැ
  ඔ
  නො
  දෙ

=head3 Wide (stand-alone)

  ජනවාර
  පෙබරවාර
  මාර්ත
  අප්‍රේල්
  මැයි
  ජූන
  ජූලි
  අගෝස්තු
  සැප්තැම්බර්
  ඔක්තෝබර්
  නොවැම්බර්
  දෙසැම්බර්

=head3 Abbreviated (stand-alone)

  ජන
  පෙබ
  මාර්ත
  අප්‍රේල
  මැය
  ජූන
  ජූල
  අගෝ
  සැප
  ඔක
  නොවැ
  දෙසැ

=head3 Narrow (stand-alone)

  ජ
  පෙ
  මා
  අ
  මැ
  ජූ
  ජූ
  අ
  සැ
  ඔ
  නො
  දෙ

=head2 Quarters

=head3 Wide (format)

  1 වන කාර්තුව
  2 වන කාර්තුව
  3 වන කාර්තුව
  4 වන කාර්තුව

=head3 Abbreviated (format)

  කාර්:1
  කාර්:2
  කාර්:3
  කාර්:4

=head3 Narrow (format)

  1
  2
  3
  4

=head3 Wide (stand-alone)

  1 වන කාර්තුව
  2 වන කාර්තුව
  3 වන කාර්තුව
  4 වන කාර්තුව

=head3 Abbreviated (stand-alone)

  කාර්:1
  කාර්:2
  කාර්:3
  කාර්:4

=head3 Narrow (stand-alone)

  1
  2
  3
  4

=head2 Eras

=head3 Wide

  ක්‍රිස්තු පූර්‍ව
  ක්‍රිස්තු වර්‍ෂ

=head3 Abbreviated

  ක්‍රි.පූ.
  ක්‍රි.ව.

=head3 Narrow

  ක්‍රි.පූ.
  ක්‍රි.ව.

=head2 Date Formats

=head3 Full

   2008-02-05T18:30:30 = අඟහරුවාදා, 2008 පෙබරවාර 5
   1995-12-22T09:05:02 = සිකුරාදා, 1995 දෙසැම්බර් 22
  -0010-09-15T04:44:23 = සෙනසුරාදා, -10 සැප්තැම්බර් 15

=head3 Long

   2008-02-05T18:30:30 = 2008 පෙබරවාර 5
   1995-12-22T09:05:02 = 1995 දෙසැම්බර් 22
  -0010-09-15T04:44:23 = -10 සැප්තැම්බර් 15

=head3 Medium

   2008-02-05T18:30:30 = 2008 පෙබ 5
   1995-12-22T09:05:02 = 1995 දෙසැ 22
  -0010-09-15T04:44:23 = -10 සැප 15

=head3 Short

   2008-02-05T18:30:30 = 2008/02/05
   1995-12-22T09:05:02 = 1995/12/22
  -0010-09-15T04:44:23 = -010/09/15

=head3 Default

   2008-02-05T18:30:30 = 2008 පෙබ 5
   1995-12-22T09:05:02 = 1995 දෙසැ 22
  -0010-09-15T04:44:23 = -10 සැප 15

=head2 Time Formats

=head3 Full

   2008-02-05T18:30:30 = 6:30:30 ප.ව. UTC
   1995-12-22T09:05:02 = 9:05:02 පෙ.ව. UTC
  -0010-09-15T04:44:23 = 4:44:23 පෙ.ව. UTC

=head3 Long

   2008-02-05T18:30:30 = 6:30:30 ප.ව. UTC
   1995-12-22T09:05:02 = 9:05:02 පෙ.ව. UTC
  -0010-09-15T04:44:23 = 4:44:23 පෙ.ව. UTC

=head3 Medium

   2008-02-05T18:30:30 = 6:30:30 ප.ව.
   1995-12-22T09:05:02 = 9:05:02 පෙ.ව.
  -0010-09-15T04:44:23 = 4:44:23 පෙ.ව.

=head3 Short

   2008-02-05T18:30:30 = 6:30 ප.ව.
   1995-12-22T09:05:02 = 9:05 පෙ.ව.
  -0010-09-15T04:44:23 = 4:44 පෙ.ව.

=head3 Default

   2008-02-05T18:30:30 = 6:30:30 ප.ව.
   1995-12-22T09:05:02 = 9:05:02 පෙ.ව.
  -0010-09-15T04:44:23 = 4:44:23 පෙ.ව.

=head2 Datetime Formats

=head3 Full

   2008-02-05T18:30:30 = අඟහරුවාදා, 2008 පෙබරවාර 5 6:30:30 ප.ව. UTC
   1995-12-22T09:05:02 = සිකුරාදා, 1995 දෙසැම්බර් 22 9:05:02 පෙ.ව. UTC
  -0010-09-15T04:44:23 = සෙනසුරාදා, -10 සැප්තැම්බර් 15 4:44:23 පෙ.ව. UTC

=head3 Long

   2008-02-05T18:30:30 = 2008 පෙබරවාර 5 6:30:30 ප.ව. UTC
   1995-12-22T09:05:02 = 1995 දෙසැම්බර් 22 9:05:02 පෙ.ව. UTC
  -0010-09-15T04:44:23 = -10 සැප්තැම්බර් 15 4:44:23 පෙ.ව. UTC

=head3 Medium

   2008-02-05T18:30:30 = 2008 පෙබ 5 6:30:30 ප.ව.
   1995-12-22T09:05:02 = 1995 දෙසැ 22 9:05:02 පෙ.ව.
  -0010-09-15T04:44:23 = -10 සැප 15 4:44:23 පෙ.ව.

=head3 Short

   2008-02-05T18:30:30 = 2008/02/05 6:30 ප.ව.
   1995-12-22T09:05:02 = 1995/12/22 9:05 පෙ.ව.
  -0010-09-15T04:44:23 = -010/09/15 4:44 පෙ.ව.

=head3 Default

   2008-02-05T18:30:30 = 2008 පෙබ 5 6:30:30 ප.ව.
   1995-12-22T09:05:02 = 1995 දෙසැ 22 9:05:02 පෙ.ව.
  -0010-09-15T04:44:23 = -10 සැප 15 4:44:23 පෙ.ව.

=head2 Available Formats

=head3 d (d)

   2008-02-05T18:30:30 = 5
   1995-12-22T09:05:02 = 22
  -0010-09-15T04:44:23 = 15

=head3 EEEd (d EEE)

   2008-02-05T18:30:30 = 5 අඟ
   1995-12-22T09:05:02 = 22 සිකු
  -0010-09-15T04:44:23 = 15 සෙන

=head3 Hm (H:mm)

   2008-02-05T18:30:30 = 18:30
   1995-12-22T09:05:02 = 9:05
  -0010-09-15T04:44:23 = 4:44

=head3 hm (h:mm a)

   2008-02-05T18:30:30 = 6:30 ප.ව.
   1995-12-22T09:05:02 = 9:05 පෙ.ව.
  -0010-09-15T04:44:23 = 4:44 පෙ.ව.

=head3 Hms (H:mm:ss)

   2008-02-05T18:30:30 = 18:30:30
   1995-12-22T09:05:02 = 9:05:02
  -0010-09-15T04:44:23 = 4:44:23

=head3 hms (h:mm:ss a)

   2008-02-05T18:30:30 = 6:30:30 ප.ව.
   1995-12-22T09:05:02 = 9:05:02 පෙ.ව.
  -0010-09-15T04:44:23 = 4:44:23 පෙ.ව.

=head3 M (L)

   2008-02-05T18:30:30 = 2
   1995-12-22T09:05:02 = 12
  -0010-09-15T04:44:23 = 9

=head3 Md (M-d)

   2008-02-05T18:30:30 = 2-5
   1995-12-22T09:05:02 = 12-22
  -0010-09-15T04:44:23 = 9-15

=head3 MEd (E, M-d)

   2008-02-05T18:30:30 = අඟ, 2-5
   1995-12-22T09:05:02 = සිකු, 12-22
  -0010-09-15T04:44:23 = සෙන, 9-15

=head3 MMM (LLL)

   2008-02-05T18:30:30 = පෙබ
   1995-12-22T09:05:02 = දෙසැ
  -0010-09-15T04:44:23 = සැප

=head3 MMMd (MMM d)

   2008-02-05T18:30:30 = පෙබ 5
   1995-12-22T09:05:02 = දෙසැ 22
  -0010-09-15T04:44:23 = සැප 15

=head3 MMMEd (E MMM d)

   2008-02-05T18:30:30 = අඟ පෙබ 5
   1995-12-22T09:05:02 = සිකු දෙසැ 22
  -0010-09-15T04:44:23 = සෙන සැප 15

=head3 MMMMd (MMMM d)

   2008-02-05T18:30:30 = පෙබරවාර 5
   1995-12-22T09:05:02 = දෙසැම්බර් 22
  -0010-09-15T04:44:23 = සැප්තැම්බර් 15

=head3 MMMMEd (E MMMM d)

   2008-02-05T18:30:30 = අඟ පෙබරවාර 5
   1995-12-22T09:05:02 = සිකු දෙසැම්බර් 22
  -0010-09-15T04:44:23 = සෙන සැප්තැම්බර් 15

=head3 ms (mm:ss)

   2008-02-05T18:30:30 = 30:30
   1995-12-22T09:05:02 = 05:02
  -0010-09-15T04:44:23 = 44:23

=head3 y (y)

   2008-02-05T18:30:30 = 2008
   1995-12-22T09:05:02 = 1995
  -0010-09-15T04:44:23 = -10

=head3 yM (y-M)

   2008-02-05T18:30:30 = 2008-2
   1995-12-22T09:05:02 = 1995-12
  -0010-09-15T04:44:23 = -10-9

=head3 yMEd (EEE, y-M-d)

   2008-02-05T18:30:30 = අඟ, 2008-2-5
   1995-12-22T09:05:02 = සිකු, 1995-12-22
  -0010-09-15T04:44:23 = සෙන, -10-9-15

=head3 yMMM (y MMM)

   2008-02-05T18:30:30 = 2008 පෙබ
   1995-12-22T09:05:02 = 1995 දෙසැ
  -0010-09-15T04:44:23 = -10 සැප

=head3 yMMMEd (EEE, y MMM d)

   2008-02-05T18:30:30 = අඟ, 2008 පෙබ 5
   1995-12-22T09:05:02 = සිකු, 1995 දෙසැ 22
  -0010-09-15T04:44:23 = සෙන, -10 සැප 15

=head3 yMMMM (y MMMM)

   2008-02-05T18:30:30 = 2008 පෙබරවාර
   1995-12-22T09:05:02 = 1995 දෙසැම්බර්
  -0010-09-15T04:44:23 = -10 සැප්තැම්බර්

=head3 yQ (y Q)

   2008-02-05T18:30:30 = 2008 1
   1995-12-22T09:05:02 = 1995 4
  -0010-09-15T04:44:23 = -10 3

=head3 yQQQ (y QQQ)

   2008-02-05T18:30:30 = 2008 කාර්:1
   1995-12-22T09:05:02 = 1995 කාර්:4
  -0010-09-15T04:44:23 = -10 කාර්:3

=head2 Miscellaneous

=head3 Prefers 24 hour time?

No

=head3 Local first day of the week

සඳුදා


=head1 SUPPORT

See L<DateTime::Locale>.

=head1 AUTHOR

Dave Rolsky <autarch@urth.org>

=head1 COPYRIGHT

Copyright (c) 2008 David Rolsky. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same
terms as Perl itself.

This module was generated from data provided by the CLDR project, see
the LICENSE.cldr in this distribution for details on the CLDR data's
license.

=cut
