# -*- coding: utf-8 -*-
#
#== Fortran90 code generator as macro
#
#Authors::   Yasuhiro MORIKAWA
#Version::   $Id: lib-rb2f90-macro.rb,v 1.1 2009-03-25 08:17:36 morikawa Exp $
#Tag Name::  $Name:  $
#Copyright:: Copyright (C) GFD Dennou Club, 2005-2009. All rights reserved.
#License::   See link:../COPYRIGHT
#
#These functions are used to generate f90 code files from Ruby code 
#files. They are expected to help ruby code approximates f90 code for
#as long as possible.
#
#[JAPANESE]
#
#これらの関数は Ruby で記述されたファイルから F90 ファイルを生成
#するための関数です. これらの関数により, できるだけ F90 コードに
#近い形で Ruby のコードを記述できることが期待されます.
#

require("lib-rb2f90-macro-intrinsic_types")
#
#== Header Comment
#
def rb2f90_header_comment
  mess = <<-"__EOF__"
! *** Caution!! ***
! 
! This file is generated from "#{File.basename($0.to_s)}" by Ruby #{RUBY_VERSION}.
! Please do not edit this file directly.
!
! [JAPANESE]
!
! ※※※ 注意!!! ※※※
!
! このファイルは "#{File.basename($0.to_s)}" から Ruby #{RUBY_VERSION}
! によって自動生成されたファイルです.
! このファイルを直接編集しませんようお願い致します.
!
__EOF__

  return mess
end

#
#== Local variable section in emacs
#
def rb2f90_emacs_readonly
  mess = "!Local Variables:\n"
  mess << <<-"__EOF__"
!mode: f90
!buffer-read-only: t
!End:
__EOF__
  return mess
end

#
#== An imitation of the ifelse function of m4.
#
#=== Usage
#
#   ifelse(string-1, string-2, equal, [not-equal])
#   ifelse(string-1, string-2, equal, [string-3, string-4, equal, ... [not-equal]])
#
# [JA]
#
#== m4 の ifelse 関数を模した関数
#
#=== 使い方
#
#    ifelse(string-1, string-2, equal, [not-equal])
#    ifelse(string-1, string-2, equal, [string-3, string-4, equal, ... [not-equal]])
# string-1 と string-2 が等しい場合, equal が返る. 等しくない場合には
# not-equal が返る. not-equal が指定されない場合には何も返らない.
#
# 後ろにさらに条件分岐を付け加えることも可能である.
# 3 つ目の例では, string-3 と 4 が指定されている.
# string-1 と 2 が異なる場合, 次に string-3 と string-4 が比較され,
# 等しい場合にはその後ろの equal が返る. 等しくない場合には
# さらに後ろの引数に制御を渡す.
#
def ifelse(*all)
  entire = Array.new
  count = -1
  one_set = Array.new(3)
  body = ""
  # 正規化 (?)
  # 3 の倍数 + 2 の場合, 最後の1つは捨てる
  if all.size.modulo(3) == 2 then
    all.pop
  end
  # 3 の倍数 + 1 の場合, 最後の1つは別途持っておく
  lastitem = Array.new
  if all.size.modulo(3) == 1 then
    lastitem << all.pop
  end
  # データの整理
  all.each{ |item|
    count += 1
    one_set[count.modulo(3)] = item
    if count.modulo(3) == 2 then
      entire << one_set.clone
      one_set.clear
      next
    end
  }
  if !lastitem.empty? then
    entire << lastitem
  end
  entire.each{ |set|
    if set.size == 1 then
      body << set[0].sub(/^\n/, '').chomp
      break
    end
    if set[0] == set[1] then
      body << set[2].sub(/^\n/, '').chomp
      break
    end
  }
  return body.chomp
end

#
#== An imitation of the forloop function in m4-doc.
#
#=== Usage
#
#   forloop(string-1, first, last, string-2)
#
# [JA]
#
#== m4-doc で紹介される forloop 関数を模した関数
#
#=== 使い方
#
#   forloop(str, first, last, body)
#
# first, last には整数を代入する.
# body を last - first 回文だけ反復して返す. 
# その際, body  内で str と同じ文字列を数値に置換する.
# 置換される数値は first から last へと 1 つづつ増加する数値である.
#
def forloop(str, first, last, body)
  rbody = ""
  repeated = nil
  for i in first..last
    rbody << "\n" if repeated
    rbody << body.sub(/^\n/, '').gsub(/#{str}/, i.to_s).chomp
    repeated = true
  end
  return rbody.chomp
end

#
#== Function as "for" for macro
#
#=== Usage
#
#   foreach(string-1, first, last, string-2)
#
# [JA]
#
#== マクロ的に利用できる for 関数
#
#=== 使い方
#
#   forloop(string-1, first, last, string-2)
#
# first, last には整数を代入する.
# string-2 を last - first 回文だけ反復して返す. 
# その際, string-2 内で string-1 と同じ文字列を数値に置換する.
# 置換される数値は first から last へと 1 つづつ増加する数値である.
#
def foreach(str, *words)
  return if words.size < 2
  body = "#{words[words.size - 1]}"
  words.pop
  rbody = ""
  repeated = nil
  words.each{ |word|
    rbody << "\n" if repeated
    rbody << body.sub(/^\n/, '').gsub(/#{str}/, word.to_s).chomp
    repeated = true
  }
  return rbody.chomp
end

#
#== Output (:), or (:,:) or ... in array definition
#
#=== Usage
#
#   array_colon(num)
#
# [JA]
#
#== 配列定義の際の (:), (:,:), ... を出力する関数
#
#=== 使い方
#
#   array_colon(num)
#
# num に正の整数を代入する. 0 ならば何も出力せず, 1 ならば (:), 
# 2 ならば (:,:), ... と出力する.
#
def array_colon(num)
  return unless num
  int = num.to_i
  return if int < 1
  return "(:)" if int < 2
  int -= 1
  body = "(:"
  int.times{
    body << ",:"
  }
  body << ")"
  return body
end

